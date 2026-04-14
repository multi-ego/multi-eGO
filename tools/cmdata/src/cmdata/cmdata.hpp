#ifndef _CMDATA_CMDATA_HPP
#define _CMDATA_CMDATA_HPP

// gromacs includes
#include <gromacs/trajectoryanalysis/topologyinformation.h>
#ifdef GMXVGE2026
#include "gromacs/utility/vec.h"
#else
#include <gromacs/math/vec.h>
#endif
#include <gromacs/pbcutil/pbc.h>
#include <gromacs/fileio/tpxio.h>

// cmdata includes
#include "io.hpp"
#include "indexing.hpp"
#include "density.hpp"
#include "mindist.hpp"
#include "xtc_frame.hpp"
#include "function_types.hpp"

// standard library imports
#include <iostream>
#include <cmath>
#include <memory>
#include <string>
#include <algorithm>
#include <functional>
#include <numeric>
#include <sstream>
#include <fstream>

// xdrfile includes
#include <xdrfile.h>
#include <xdrfile_xtc.h>

namespace cmdata
{

class CMData
{
private:
  // general fields
  float cutoff_, mol_cutoff_, mcut2_, cut_sig_2_;
  int nskip_, n_x_;
  float dt_, t_begin_, t_end_;
  gmx_mtop_t *mtop_;
  rvec *xcm_ = nullptr;

  // molecule number fields
  int nindex_;
  gmx::RangePartitioning mols_;
  std::vector<int> natmol2_;
  std::vector<int> mol_id_;
  std::vector<int> num_mol_unique_;
  std::vector<float> inv_num_mol_unique_;
  std::vector<float> inv_num_mol_;
  std::string bkbn_H_;

  // atom_active_[mt][a] == true when atom a of molecule type mt should be
  // included in distance calculations (i.e. not a skipped hydrogen).
  // Precomputed once in initAnalysis() to avoid per-frame topology lookups.
  std::vector<std::vector<bool>> atom_active_;

  // weights fields
  float weights_sum_;
  std::string weights_path_;
  std::vector<float> weights_;

  // pbc fields
  bool no_pbc_;
  t_pbc *pbc_;
  PbcType pbcType_;

  // frame fields
  cmdata::xtc::Frame *frame_;
  XDRFILE *trj_;

  // density fields
  float dx_;
  std::size_t n_bins_;
  std::vector<float> density_bins_;
  std::vector<std::vector<int>> cross_index_;

  using cmdata_matrix = std::vector<std::vector<std::vector<std::vector<float>>>>;
  cmdata_matrix interm_same_mat_density_;
  cmdata_matrix interm_cross_mat_density_;
  cmdata_matrix intram_mat_density_;
  cmdata_matrix interm_same_maxcdf_mol_;
  cmdata_matrix interm_cross_maxcdf_mol_;

  // temporary containers for maxcdf operations
  std::vector<std::vector<float>> frame_same_mat_;
  std::vector<std::vector<float>> frame_cross_mat_;

  // mode selection, booleans and functions
  std::string mode_;
  bool intra_ = false, same_ = false, cross_ = false;
  bool h5_ = false;

  // No std::function members for the density routines — molecule_routine is
  // templated on FIntra/FSame/FCross so the compiler can inline all three.

  template<typename FIntra, typename FSame, typename FCross>
  static void molecule_routine(
    const int i, const int nindex_, t_pbc *pbc, rvec *x, const std::vector<float> &inv_num_mol_, const float cut_sig_2_,
    const std::vector<int> &natmol2_, const std::vector<int> &num_mol_unique_, const std::vector<int> &mol_id_,
    const std::vector<std::vector<int>> &cross_index_, const std::vector<float> &density_bins_, const float mcut2_,
    rvec *xcm_, const gmx::RangePartitioning &mols_,
    const std::vector<std::vector<bool>> &atom_active_,
    std::vector<std::vector<float>> &frame_same_mat_,
    cmdata_matrix &intram_mat_density_, cmdata_matrix &interm_same_mat_density_, std::vector<std::vector<float>> &frame_cross_mat_,
    cmdata_matrix &interm_cross_mat_density_,
    const FIntra &f_intra_mol_, const FSame &f_inter_mol_same_, const FCross &f_inter_mol_cross_,
    float weight
  )
  {
    int tmp_i = 0;
    std::size_t mol_i = i, mol_j = 0;
    while ( static_cast<int>(mol_i) - num_mol_unique_[tmp_i] >= 0 )
    {
      mol_i -= num_mol_unique_[tmp_i];
      tmp_i++;
      if (tmp_i == static_cast<int>(num_mol_unique_.size())) break;
    }
    if (mol_i == static_cast<std::size_t>(num_mol_unique_[mol_id_[i]])) mol_i = 0;

    const int mt_i = mol_id_[i];
    /* Loop over molecules  */
    for (int j = 0; j < nindex_; j++)
    {
      if (j != 0)
        if (mol_j == static_cast<std::size_t>(num_mol_unique_[mol_id_[j-1]])) mol_j = 0;

      /* intermolecular interactions are evaluated only among neighbour molecules */
      if (i != j)
      {
        rvec dx;
        if (pbc != nullptr) pbc_dx(pbc, xcm_[i], xcm_[j], dx);
        else rvec_sub(xcm_[i], xcm_[j], dx);
        float dx2 = iprod(dx, dx);
        if (dx2 > mcut2_) continue;
      }
      /* for molecules of different specie we fill half a matrix */
      if (mol_id_[i] != mol_id_[j] && j < i) continue;
      std::size_t a_i = 0;

      const int mt_j = mol_id_[j];
      /* cycle over the atoms of a molecule i */
      for (std::size_t ii = mols_.block(i).begin(); ii < mols_.block(i).end(); ii++)
      {
        std::size_t a_j = 0;
        if (!atom_active_[mt_i][a_i])
        {
          a_i++;
          continue;
        }
        /* cycle over the atoms of a molecule j */
        for (std::size_t jj = mols_.block(j).begin(); jj < mols_.block(j).end(); jj++)
        {
          if (!atom_active_[mt_j][a_j])
          {
            a_j++;
            continue;
          }
          std::size_t delta  = a_i - a_j;
          rvec sym_dx;
          if (pbc != nullptr) pbc_dx(pbc, x[ii], x[jj], sym_dx);
          else rvec_sub(x[ii], x[jj], sym_dx);
          float dx2 = iprod(sym_dx, sym_dx);
          if (i==j)
          {
            if (dx2 < cut_sig_2_)
            { // intra molecule species
              f_intra_mol_(i, a_i, a_j, dx2, weight, mol_id_, natmol2_, density_bins_, inv_num_mol_, intram_mat_density_);
            }
          }
          else
          {
            if (mol_id_[i]==mol_id_[j])
            { // inter same molecule specie
              if (dx2 < cut_sig_2_)
              {
                f_inter_mol_same_(
                  i, mol_i, a_i, a_j, dx2, weight, mol_id_, natmol2_, density_bins_, frame_same_mat_, interm_same_mat_density_
                );
              }
              if (delta != 0) {
                // this is to account for inversion atom/molecule
                if (pbc != nullptr) pbc_dx(pbc, x[ii-delta], x[jj+delta], sym_dx);
                else rvec_sub(x[ii-delta], x[jj+delta], sym_dx);
                dx2 = iprod(sym_dx, sym_dx);
                if (dx2 < cut_sig_2_)
                {
                  f_inter_mol_same_(
                    i, mol_i, a_i, a_j, dx2, weight, mol_id_, natmol2_, density_bins_, frame_same_mat_, interm_same_mat_density_
                  );
                }
              }
            }
            else
            { // inter cross molecule species
              if (dx2 < cut_sig_2_)
              {
                f_inter_mol_cross_(
                  i, j, mol_i, mol_j, a_i, a_j, dx2, weight, mol_id_, natmol2_, cross_index_, density_bins_, num_mol_unique_, frame_cross_mat_, interm_cross_mat_density_
                );
              }
            }
          }
          ++a_j;
        }
        ++a_i;
      }
      ++mol_j;
    }
  }

public:
  CMData(
    const std::string &top_path, const std::string &traj_path,
    float cutoff, float mol_cutoff, int nskip,
    int dt, const std::string &mode, const std::string &bkbn_H, const std::string &weights_path,
    bool no_pbc, float t_begin, float t_end, bool h5
  ) : cutoff_(cutoff), mol_cutoff_(mol_cutoff), nskip_(nskip), dt_(dt), t_begin_(t_begin), t_end_(t_end),
      bkbn_H_(bkbn_H), weights_path_(weights_path), no_pbc_(no_pbc), mode_(mode), h5_(h5)
  {
    matrix boxtop_;
    mtop_ = (gmx_mtop_t*)malloc(sizeof(gmx_mtop_t));
    int natoms;
    pbcType_ = read_tpx(top_path.c_str(), nullptr, boxtop_, &natoms, nullptr, nullptr, mtop_);

    if (no_pbc_)
    {
      pbc_ = nullptr;
    }
    else
    {
      pbc_ = (t_pbc*)malloc(sizeof(t_pbc));
      set_pbc(pbc_, pbcType_, boxtop_);
    }

    int natom;
    long unsigned int nframe;
    int64_t *offsets;

    frame_ = (cmdata::xtc::Frame*)malloc(sizeof(cmdata::xtc::Frame));
    std::cout << "Reading trajectory file " << traj_path << std::endl;
    read_xtc_header(traj_path.c_str(), &natom, &nframe, &offsets);
    *frame_ = cmdata::xtc::Frame(natom);
    frame_->nframe  = nframe;
    frame_->offsets = offsets;

    trj_ = xdrfile_open(traj_path.c_str(), "r");
    initAnalysis();
  }

  ~CMData()
  {
    free(frame_->x);
    free(frame_->offsets);
    free(frame_);
    free(mtop_);
    if (xcm_ != nullptr) free(xcm_);
  }

  void initAnalysis()
  /**
   * @brief Initializes the analysis by setting up the molecule partitioning and the mode selection
   *
   * @todo Check if molecule block is empty
  */
  {
    n_x_ = 0;

    // get the number of atoms per molecule
    // equivalent to mols_ = gmx:gmx_mtop_molecules(*top.mtop());
    for (const gmx_molblock_t &molb : mtop_->molblock)
    {
      int natm_per_mol = mtop_->moltype[molb.type].atoms.nr;
      for (int i = 0; i < molb.nmol; i++) mols_.appendBlock(natm_per_mol);
      num_mol_unique_.push_back(molb.nmol);
    }
    // number of molecules
    nindex_ = mols_.numBlocks();

    printf("Evaluating mode selection:\n");
    std::string tmp_mode;
    std::stringstream modestream{ mode_ };
    while (std::getline(modestream, tmp_mode, '+'))
    {
      if (tmp_mode == std::string("intra"))
      {
        intra_ = true;
      }
      else if (tmp_mode == std::string("same"))
      {
        same_ = true;
      }
      else if (tmp_mode == std::string("cross"))
      {
        cross_ = true;
      }
      else
      {
        printf("Wrong mode: %s\nMode must be one from: intra, same, cross. Use + to concatenate more than one, i.e. intra+cross\n", tmp_mode.c_str());
        exit(1);
      }
      printf(" - found %s\n", tmp_mode.c_str());
    }

    int mol_id = 0;
    int molb_index = 0;
    for ( auto i : num_mol_unique_ )
    {
      natmol2_.push_back(mols_.block(molb_index).end() - mols_.block(molb_index).begin());
      inv_num_mol_unique_.push_back(1. / static_cast<float>(i));
      for ( int j = 0; j < i; j++ )
      {
        mol_id_.push_back(mol_id);
        inv_num_mol_.push_back(1. / static_cast<float>(i));
      }
      mol_id++;
      molb_index += i;
    }

    // Precompute per-molecule-type atom activity mask.
    // Walk the first molecule of each type and record which atoms to include.
    {
      int molb = 0;
      int mol_first = 0; // global index of the first molecule of each type
      atom_active_.resize(natmol2_.size());
      for (std::size_t mt = 0; mt < natmol2_.size(); mt++)
      {
        atom_active_[mt].resize(natmol2_[mt]);
        const char *atomname;
        for (int a = 0; a < natmol2_[mt]; a++)
        {
          int global_atom = mols_.block(mol_first).begin() + a;
          mtopGetAtomAndResidueName(*mtop_, global_atom, &molb, &atomname, nullptr, nullptr, nullptr);
          // Skip non-backbone hydrogens.  "H" (AMBER/GROMACS) and "HN" (CHARMM)
          // are always treated as backbone H; --bkbn_H adds a third custom name.
          {
            const std::string aname(atomname);
            const bool is_H = (aname[0] == 'H');
            const bool is_bkbn = (aname == "H" || aname == "HN" ||
                                  (!bkbn_H_.empty() && aname == bkbn_H_));
            atom_active_[mt][a] = !(is_H && !is_bkbn);
          }
        }
        mol_first += num_mol_unique_[mt]; // advance to first mol of next type
      }
    }

    printf("Number of different molecules %lu\n", natmol2_.size());
    bool check_same = false;
    for(std::size_t i=0; i<natmol2_.size();i++) {
      printf("mol %lu num %u size %u\n", i, num_mol_unique_[i], natmol2_[i]);
      if(num_mol_unique_[i]>1) check_same = true;
    }
    if(!check_same && same_) same_ = false;
    if(natmol2_.size()<2) cross_ = false;
    if(nindex_>1 && (same_ || cross_))
    {
      printf("\n\n::::::::::::WARNING::::::::::::\nMore than 1 molecule found in the system.\nFix pbc before running cmdata using pbc mol\n");
      printf(":::::::::::::::::::::::::::::::\n\n");
    }
    if (same_)
    {
      std::cout << ":: activating intermat same calculations" << std::endl;
      interm_same_mat_density_.resize(natmol2_.size());
      interm_same_maxcdf_mol_.resize(natmol2_.size());
    }
    if (cross_)
    {
      std::cout << ":: activating intermat cross calculations" << std::endl;
      interm_cross_mat_density_.resize((natmol2_.size() * (natmol2_.size() - 1)) / 2);
      interm_cross_maxcdf_mol_.resize((natmol2_.size() * (natmol2_.size() - 1)) / 2);
    }
    if (intra_)
    {
      std::cout << " :: activating intramat calculations" << std::endl;
      intram_mat_density_.resize(natmol2_.size());
    }

    density_bins_.resize(cmdata::indexing::n_bins(cutoff_));
    for (std::size_t i = 0; i < density_bins_.size(); i++)
      density_bins_[i] = cutoff_ / static_cast<float>(density_bins_.size()) * static_cast<float>(i) + cutoff_ / static_cast<float>(density_bins_.size() * 2);

    int cross_count = 0;
    if (cross_) cross_index_.resize(natmol2_.size(), std::vector<int>(natmol2_.size(), 0));
    for ( std::size_t i = 0; i < natmol2_.size(); i++ )
    {
      if (same_)
      {
        interm_same_mat_density_[i].resize(natmol2_[i], std::vector<std::vector<float>>(natmol2_[i], std::vector<float>(cmdata::indexing::n_bins(cutoff_), 0)));
        interm_same_maxcdf_mol_[i].resize(natmol2_[i], std::vector<std::vector<float>>(natmol2_[i], std::vector<float>(cmdata::indexing::n_bins(cutoff_), 0)));
      }
      if (intra_) intram_mat_density_[i].resize(natmol2_[i], std::vector<std::vector<float>>(natmol2_[i], std::vector<float>(cmdata::indexing::n_bins(cutoff_), 0)));
      for ( std::size_t j = i + 1; j < natmol2_.size() && cross_; j++ )
      {
        interm_cross_mat_density_[cross_count].resize(natmol2_[i], std::vector<std::vector<float>>(natmol2_[j], std::vector<float>(cmdata::indexing::n_bins(cutoff_), 0)));
        interm_cross_maxcdf_mol_[cross_count].resize(natmol2_[i], std::vector<std::vector<float>>(natmol2_[j], std::vector<float>(cmdata::indexing::n_bins(cutoff_), 0)));
        cross_index_[i][j] = cross_count;
        cross_count++;
      }
    }

    n_bins_ = cmdata::indexing::n_bins(cutoff_);
    dx_ = cutoff_ / static_cast<float>(n_bins_);

    mcut2_ = mol_cutoff_ * mol_cutoff_;
    cut_sig_2_ = (cutoff_ + 0.02) * (cutoff_ + 0.02);
    xcm_ = (rvec*)malloc(nindex_ * sizeof(rvec));

    if (same_) frame_same_mat_.resize(natmol2_.size());
    if (cross_) frame_cross_mat_.resize(cross_index_.size());
    for ( std::size_t i = 0; i < natmol2_.size(); i++ )
    {
      if (same_) frame_same_mat_[i].resize(natmol2_[i] * natmol2_[i] * num_mol_unique_[i], 0);
      for ( std::size_t j = i+1; j < natmol2_.size() && cross_; j++ )
      {
        frame_cross_mat_[cross_index_[i][j]].resize(natmol2_[i] * natmol2_[j] * num_mol_unique_[i] * num_mol_unique_[j], 0);
      }
    }

    if (weights_path_ != "")
    {
      printf("Weights file provided. Reading weights from %s\n", weights_path_.c_str());
      weights_ = cmdata::io::read_weights_file(weights_path_);
      printf("Found %li frame weights in file\n", weights_.size());
      float w_sum = std::accumulate(std::begin(weights_), std::end(weights_), 0.0, std::plus<>());
      printf("Sum of weights amounts to %lf\n", w_sum);
      weights_sum_ = 0.;
    }

    printf("Finished preprocessing.\nStarting frame-by-frame analysis.\n");
  }

  void run()
  {
    std::cout << "Running frame-by-frame analysis" << std::endl;

    // No-op callables used for disabled modes — stateless lambdas, zero overhead.
    auto no_intra = [](auto&&...) {};
    auto no_same  = [](auto&&...) {};
    auto no_cross = [](auto&&...) {};

    using I = decltype(cmdata::density::intra_mol_routine)*;
    using S = decltype(cmdata::density::inter_mol_same_routine)*;
    using C = decltype(cmdata::density::inter_mol_cross_routine)*;
    I fi = cmdata::density::intra_mol_routine;
    S fs = cmdata::density::inter_mol_same_routine;
    C fc = cmdata::density::inter_mol_cross_routine;

    // do_run instantiates molecule_routine with the correct concrete function types
    // so the compiler can inline all three paths.
    auto do_run = [this](auto f_intra, auto f_same, auto f_cross)
    {
      int frnr = 0;
      float progress = 0.0, new_progress = 0.0;
      cmdata::io::print_progress_bar(progress);
      while (frame_->read_next_frame(trj_, no_pbc_, pbcType_, pbc_) == exdrOK)
      {
        new_progress = static_cast<float>(frnr) / static_cast<float>(frame_->nframe);
        if (new_progress - progress > 0.01)
        {
          progress = new_progress;
          cmdata::io::print_progress_bar(progress);
        }
        if ((frame_->time >= t_begin_ && (t_end_ < 0 || frame_->time <= t_end_)) &&
            (dt_ == 0 || std::fmod(frame_->time, dt_) == 0) && (nskip_ == 0 || std::fmod(frnr, nskip_) == 0))
        {
          float weight = 1.0;
          if (!weights_.empty())
          {
            weight = weights_[frnr];
            weights_sum_ += weight;
          }
          // TODO: replace magic sentinel 100.f (nm) with a named constant, e.g.:
          //   static constexpr float LARGE_DIST = 100.f;
          /* resetting the per-frame minimum-distance matrix to a large sentinel value */
          for (auto &v : frame_same_mat_) std::fill(v.begin(), v.end(), 100.f);
          for (auto &v : frame_cross_mat_) std::fill(v.begin(), v.end(), 100.f);

          /* compute molecule centres of mass */
          for (int i = 0; i < nindex_; i++)
          {
            clear_rvec(xcm_[i]);
            float tm = 0.;
            for (int ii = mols_.block(i).begin(); ii < mols_.block(i).end(); ii++)
            {
              for (int m = 0; m < DIM; m++) xcm_[i][m] += frame_->x[ii][m];
              tm += 1.0;
            }
            for (int m = 0; m < DIM; m++) xcm_[i][m] /= tm;
          }

          /* per-molecule distance computation */
          for (int i = 0; i < nindex_; i++)
            molecule_routine(i, nindex_, pbc_, frame_->x, inv_num_mol_,
              cut_sig_2_, natmol2_, num_mol_unique_, mol_id_, cross_index_,
              density_bins_, mcut2_, xcm_, mols_, atom_active_,
              frame_same_mat_,
              intram_mat_density_, interm_same_mat_density_, frame_cross_mat_,
              interm_cross_mat_density_,
              f_intra, f_same, f_cross, weight);

          /* accumulate per-frame minimum distances into KDE histograms */
          if (same_)
            cmdata::mindist::mindist_same(density_bins_, num_mol_unique_, natmol2_,
              frame_same_mat_, interm_same_maxcdf_mol_, weight);
          if (cross_)
            cmdata::mindist::mindist_cross(natmol2_, cross_index_, density_bins_, num_mol_unique_,
              frame_cross_mat_, interm_cross_maxcdf_mol_, weight);

          ++n_x_;
        }
        ++frnr;
      }
      cmdata::io::print_progress_bar(1.0);
    };

    if      ( intra_ &&  same_ &&  cross_) do_run(fi,       fs,       fc      );
    else if ( intra_ &&  same_ && !cross_) do_run(fi,       fs,       no_cross);
    else if ( intra_ && !same_ &&  cross_) do_run(fi,       no_same,  fc      );
    else if ( intra_ && !same_ && !cross_) do_run(fi,       no_same,  no_cross);
    else if (!intra_ &&  same_ &&  cross_) do_run(no_intra, fs,       fc      );
    else if (!intra_ &&  same_ && !cross_) do_run(no_intra, fs,       no_cross);
    else if (!intra_ && !same_ &&  cross_) do_run(no_intra, no_same,  fc      );
    else                                   do_run(no_intra, no_same,  no_cross);
  }

  void process_data()
  {
    std::cout << "\nFinished frame-by-frame analysis\n";
    std::cout << "Analyzed " << n_x_ << " frames\n";
    std::cout << "Normalizing data... " << std::endl;
    // normalisations
    float norm = ( weights_.empty() ) ? 1. / n_x_ : 1. / weights_sum_;

    using ftype_norm = cmdata::ftypes::function_traits<decltype(&cmdata::density::normalize_histo)>;
    std::function<ftype_norm::signature> f_empty = cmdata::ftypes::do_nothing<ftype_norm>();

    std::function<ftype_norm::signature> normalize_intra = (intra_) ? cmdata::density::normalize_histo : f_empty;
    std::function<ftype_norm::signature> normalize_inter_same = (same_) ? cmdata::density::normalize_histo : f_empty;
    std::function<ftype_norm::signature> normalize_inter_cross = (cross_) ? cmdata::density::normalize_histo : f_empty;

    for (std::size_t i = 0; i < natmol2_.size(); i++)
    {
      for (int ii = 0; ii < natmol2_[i]; ii++)
      {
        for (int jj = ii; jj < natmol2_[i]; jj++)
        {
          float inv_num_mol_same = inv_num_mol_unique_[i];
          normalize_inter_same(i, ii, jj, norm, inv_num_mol_same, interm_same_maxcdf_mol_);
          normalize_inter_same(i, ii, jj, norm, 1.0, interm_same_mat_density_);
          normalize_intra(i, ii, jj, norm, 1.0, intram_mat_density_);

          float sum = 0.0;
          for ( std::size_t k = (same_) ? 0 : cmdata::indexing::n_bins(cutoff_); k < cmdata::indexing::n_bins(cutoff_); k++ )
          {
            sum+= dx_ * interm_same_maxcdf_mol_[i][ii][jj][k];
            if (sum > 1.0) sum=1.0;
            interm_same_maxcdf_mol_[i][ii][jj][k] = sum;
          }
          if (same_) interm_same_mat_density_[i][jj][ii] = interm_same_mat_density_[i][ii][jj];
          if (same_) interm_same_maxcdf_mol_[i][jj][ii] = interm_same_maxcdf_mol_[i][ii][jj];
          if (intra_) intram_mat_density_[i][jj][ii] = intram_mat_density_[i][ii][jj];
        }
      }
      for (std::size_t j = i + 1; j < natmol2_.size() && cross_; j++)
      {
        for (int ii = 0; ii < natmol2_[i]; ii++)
        {
          for (int jj = 0; jj < natmol2_[j]; jj++)
          {
            float inv_num_mol_cross = inv_num_mol_unique_[i];
            normalize_inter_cross(cross_index_[i][j], ii, jj, norm, 1.0, interm_cross_mat_density_);
            normalize_inter_cross(cross_index_[i][j], ii, jj, norm, inv_num_mol_cross, interm_cross_maxcdf_mol_);

            float sum = 0.0;
            for ( std::size_t k = (cross_) ? 0 : cmdata::indexing::n_bins(cutoff_); k < cmdata::indexing::n_bins(cutoff_); k++ )
            {
              sum += dx_ * interm_cross_maxcdf_mol_[cross_index_[i][j]][ii][jj][k];
              if (sum > 1.0) sum = 1.0;
              interm_cross_maxcdf_mol_[cross_index_[i][j]][ii][jj][k] = sum;
            }
          }
        }
      }
    }
  }

  void write_output( const std::string &output_prefix )
  {
    std::cout << "Writing data... " << std::endl;
    using ftype_write_intra = cmdata::ftypes::function_traits<decltype(&cmdata::io::f_write_intra)>;
    using ftype_write_inter_same = cmdata::ftypes::function_traits<decltype(&cmdata::io::f_write_inter_same)>;
    using ftype_write_inter_cross = cmdata::ftypes::function_traits<decltype(&cmdata::io::f_write_inter_cross)>;
    std::function<ftype_write_intra::signature> write_intra = cmdata::ftypes::do_nothing<ftype_write_intra>();
    std::function<ftype_write_inter_same::signature> write_inter_same = cmdata::ftypes::do_nothing<ftype_write_inter_same>();
    std::function<ftype_write_inter_cross::signature> write_inter_cross = cmdata::ftypes::do_nothing<ftype_write_inter_cross>();

    if (intra_)
    {
      #ifdef USE_HDF5
      if(h5_) write_intra = cmdata::io::f_write_intra_HDF5;
      else write_intra = cmdata::io::f_write_intra;
      #else
      write_intra = cmdata::io::f_write_intra;
      #endif
    }
    if (same_)
    {
      #ifdef USE_HDF5
      if(h5_) write_inter_same = cmdata::io::f_write_inter_same_HDF5;
      else write_inter_same = cmdata::io::f_write_inter_same;
      #else
      write_inter_same = cmdata::io::f_write_inter_same;
      #endif
    }
    if (cross_)
    {
      #ifdef USE_HDF5
      if(h5_) write_inter_cross = cmdata::io::f_write_inter_cross_HDF5;
      else write_inter_cross = cmdata::io::f_write_inter_cross;
      #else
      write_inter_cross = cmdata::io::f_write_inter_cross;
      #endif
    }

    for (std::size_t i = 0; i < natmol2_.size(); i++)
    {
      std::cout << "Writing data for molecule " << i << "..." << std::endl;
      cmdata::io::print_progress_bar(0.0);
      float progress = 0.0, new_progress = 0.0;

      for (int ii = 0; ii < natmol2_[i]; ii++)
      {
        new_progress = static_cast<float>(ii) / static_cast<float>(natmol2_[i]);
        if (new_progress - progress > 0.01)
        {
          progress = new_progress;
          cmdata::io::print_progress_bar(progress);
        }
        write_intra(output_prefix, i, ii, density_bins_, natmol2_, intram_mat_density_);
        write_inter_same(output_prefix, i, ii, density_bins_, natmol2_, interm_same_mat_density_, interm_same_maxcdf_mol_);
      }
      if(cross_)
      {
        for (std::size_t j = i + 1; j < natmol2_.size(); j++)
        {
          for (int ii = 0; ii < natmol2_[i]; ii++)
          {
            write_inter_cross(output_prefix, i, j, ii, density_bins_, natmol2_, cross_index_, interm_cross_mat_density_, interm_cross_maxcdf_mol_);
          }
        }
      }
    }

    cmdata::io::print_progress_bar(1.0);
    std::cout << "\nFinished!" << std::endl;
  }
};

} // namespace cmdata

#endif // _CMDATA_CMDATA_HPP
