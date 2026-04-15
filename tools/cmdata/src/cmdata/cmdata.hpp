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
#include "frame.hpp"
#include "molfile_support.hpp"

// standard library imports
#include <iostream>
#include <cmath>
#include <string>
#include <algorithm>
#include <numeric>
#include <sstream>


namespace cmdata
{

class CMData
{
private:
  // general fields
  static constexpr float LARGE_DIST = 100.f; // sentinel for per-frame minimum distance reset
  float cutoff_, mcut2_, cut_sig_2_;
  int nskip_, n_x_;
  float dt_, t_begin_, t_end_;
  // mtop_ is non-null only for TPR topology; null for structure-file topology.
  gmx_mtop_t *mtop_ = nullptr;
  // atom_names_ holds names read from non-TPR structure files (parallel to global atom index).
  std::vector<std::string> atom_names_;
  // nontpr_chain_sizes_: atoms per chain, in order, parsed from the structure file.
  // Consecutive chains of the same size are treated as one molecule type (nmol > 1).
  std::vector<int> nontpr_chain_sizes_;
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
  float weights_sum_ = 0.f;
  std::string weights_path_;
  std::vector<float> weights_;

  // pbc fields
  bool no_pbc_;
  t_pbc *pbc_;
  PbcType pbcType_;

  // frame fields
  cmdata::traj::Frame *frame_;

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

  // temporary containers for per-frame minimum distances
  std::vector<std::vector<float>> frame_same_mat_;
  std::vector<std::vector<float>> frame_cross_mat_;

  // mode selection
  std::string mode_;
  bool intra_ = false, same_ = false, cross_ = false;
  bool h5_ = false;

  // molecule_routine is templated on FIntra/FSame/FCross so the compiler can
  // inline all three density paths and the no-op lambdas for disabled modes.
  template<typename FIntra, typename FSame, typename FCross>
  void molecule_routine(
    int i, t_pbc *pbc, rvec *x,
    const FIntra &f_intra_mol_, const FSame &f_inter_mol_same_, const FCross &f_inter_mol_cross_,
    float weight
  )
  {
    std::size_t tmp_i = 0;
    std::size_t mol_i = i, mol_j = 0;
    while ( mol_i >= static_cast<std::size_t>(num_mol_unique_[tmp_i]) )
    {
      mol_i -= num_mol_unique_[tmp_i];
      tmp_i++;
      if (tmp_i == num_mol_unique_.size()) break;
    }
    if (mol_i == static_cast<std::size_t>(num_mol_unique_[mol_id_[i]])) mol_i = 0;

    const int mt_i = mol_id_[i];
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
        if (iprod(dx, dx) > mcut2_) continue;
      }
      /* for molecules of different specie we fill half a matrix */
      if (mol_id_[i] != mol_id_[j] && j < i) continue;

      const int mt_j = mol_id_[j];
      const std::size_t ii_begin = static_cast<std::size_t>(mols_.block(i).begin());
      const std::size_t ii_end   = static_cast<std::size_t>(mols_.block(i).end());
      const std::size_t jj_begin = static_cast<std::size_t>(mols_.block(j).begin());
      const std::size_t jj_end   = static_cast<std::size_t>(mols_.block(j).end());
      std::size_t a_i = 0;
      for (std::size_t ii = ii_begin; ii < ii_end; ii++)
      {
        std::size_t a_j = 0;
        if (!atom_active_[mt_i][a_i]) { a_i++; continue; }
        for (std::size_t jj = jj_begin; jj < jj_end; jj++)
        {
          if (!atom_active_[mt_j][a_j]) { a_j++; continue; }
          rvec sym_dx;
          if (pbc != nullptr) pbc_dx(pbc, x[ii], x[jj], sym_dx);
          else rvec_sub(x[ii], x[jj], sym_dx);
          float dx2 = iprod(sym_dx, sym_dx);
          if (i == j)
          {
            if (dx2 < cut_sig_2_)
              f_intra_mol_(i, a_i, a_j, dx2, weight, mol_id_, natmol2_, density_bins_, inv_num_mol_, intram_mat_density_);
          }
          else if (mol_id_[i] == mol_id_[j])
          {
            if (dx2 < cut_sig_2_)
              f_inter_mol_same_(i, mol_i, a_i, a_j, dx2, weight, mol_id_, natmol2_, density_bins_, frame_same_mat_, interm_same_mat_density_);
            if (a_i != a_j)
            {
              // Account for atom/molecule index inversion (symmetric contribution).
              // Use explicit index arithmetic to avoid unsigned underflow when a_i < a_j.
              // ii = ii_begin + a_i and jj = jj_begin + a_j always hold.
              if (pbc != nullptr) pbc_dx(pbc, x[ii_begin + a_j], x[jj_begin + a_i], sym_dx);
              else rvec_sub(x[ii_begin + a_j], x[jj_begin + a_i], sym_dx);
              dx2 = iprod(sym_dx, sym_dx);
              if (dx2 < cut_sig_2_)
                f_inter_mol_same_(i, mol_i, a_i, a_j, dx2, weight, mol_id_, natmol2_, density_bins_, frame_same_mat_, interm_same_mat_density_);
            }
          }
          else
          {
            if (dx2 < cut_sig_2_)
              f_inter_mol_cross_(i, j, mol_i, mol_j, a_i, a_j, dx2, weight, mol_id_, natmol2_, cross_index_, density_bins_, num_mol_unique_, frame_cross_mat_, interm_cross_mat_density_);
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
  ) : cutoff_(cutoff), mcut2_(mol_cutoff * mol_cutoff), nskip_(nskip), dt_(dt),
      t_begin_(t_begin), t_end_(t_end), bkbn_H_(bkbn_H), weights_path_(weights_path),
      no_pbc_(no_pbc), mode_(mode), h5_(h5)
  {
    const std::string top_ext = [&]{
      const std::size_t d = top_path.rfind('.');
      return d != std::string::npos ? top_path.substr(d + 1) : std::string{};
    }();

    if (top_ext == "tpr")
    {
      matrix boxtop_;
      mtop_ = (gmx_mtop_t*)malloc(sizeof(gmx_mtop_t));
      int natoms;
      pbcType_ = read_tpx(top_path.c_str(), nullptr, boxtop_, &natoms, nullptr, nullptr, mtop_);
      if (no_pbc_)
        pbc_ = nullptr;
      else
      {
        pbc_ = (t_pbc*)malloc(sizeof(t_pbc));
        set_pbc(pbc_, pbcType_, boxtop_);
      }
    }
    else
    {
      // Structure-file topology (pdb, gro, …): read atom names via molfile,
      // treat the whole system as one molecule, disable PBC.
      printf("Note: topology is not a .tpr — treating all atoms as one molecule, PBC disabled.\n");
      const std::string ext = top_ext.empty() ? "" : top_ext;
      molfile_plugin_t *tplugin = cmdata::get_molfile_plugin(ext);
      int natom_top = 0;
      void *thandle = tplugin->open_file_read(top_path.c_str(), ext.c_str(), &natom_top);
      if (!thandle)
        throw std::runtime_error("Cannot open topology file: " + top_path);
      atom_names_.resize(natom_top);
      if (tplugin->read_structure)
      {
        std::vector<molfile_atom_t> mf_atoms(natom_top);
        int optflags = 0;
        if (tplugin->read_structure(thandle, &optflags, mf_atoms.data()) == MOLFILE_SUCCESS)
        {
          for (int i = 0; i < natom_top; i++)
            atom_names_[i] = std::string(mf_atoms[i].name);
          // Split into chains by chain ID so each chain becomes a molecule instance.
          // Consecutive chains of the same atom count are treated as the same type.
          if (natom_top > 0)
          {
            char cur_chain = mf_atoms[0].chain[0];
            int  cur_start = 0;
            for (int i = 1; i <= natom_top; i++)
            {
              char c = (i < natom_top) ? mf_atoms[i].chain[0] : '\0';
              if (c != cur_chain)
              {
                nontpr_chain_sizes_.push_back(i - cur_start);
                cur_chain = c;
                cur_start = i;
              }
            }
          }
        }
      }
      tplugin->close_file_read(thandle);
      mtop_    = nullptr;
      pbcType_ = PbcType::No;
      pbc_     = nullptr;
      if (nontpr_chain_sizes_.empty())
        printf("Note: no chain IDs found — treating all %d atoms as one molecule.\n", natom_top);
      else
        printf("Note: detected %zu chain(s) from topology.\n", nontpr_chain_sizes_.size());
    }

    std::cout << "Reading trajectory file " << traj_path << std::endl;
    {
      const std::size_t dot = traj_path.rfind('.');
      if (dot == std::string::npos)
        throw std::runtime_error("Cannot determine trajectory format: no extension in '" + traj_path + "'");
      const std::string ext = traj_path.substr(dot + 1);
      molfile_plugin_t *plugin = cmdata::get_molfile_plugin(ext);
      int natom = 0;
      void *handle = plugin->open_file_read(traj_path.c_str(), ext.c_str(), &natom);
      if (!handle)
        throw std::runtime_error("Cannot open trajectory file: " + traj_path);
      // For non-TPR topologies the atom count comes from the structure file;
      // verify it matches what the trajectory says.
      if (!mtop_ && !atom_names_.empty() && natom != static_cast<int>(atom_names_.size()))
        throw std::runtime_error(
          "Atom count mismatch: topology has " + std::to_string(atom_names_.size()) +
          " atoms but trajectory has " + std::to_string(natom));
      // GROMACS formats (xtc/trr/gro/g96) are compiled with MOLFILE_NATIVE_NM,
      // so their coordinates are already in nm — no unit conversion needed.
      // PDB coordinates are in Angstrom and need *0.1 to reach nm.
      const bool native_nm = (ext == "xtc" || ext == "trr" || ext == "gro" || ext == "g96");
      const double coord_scale = native_nm ? 1.0 : 0.1;
      frame_ = (cmdata::traj::Frame*)malloc(sizeof(cmdata::traj::Frame));
      *frame_ = cmdata::traj::Frame(natom, plugin, handle, coord_scale);
    }
    initAnalysis();
  }

  ~CMData()
  {
    if (frame_->mf_handle && frame_->mf_plugin)
      frame_->mf_plugin->close_file_read(frame_->mf_handle);
    free(frame_->mf_coords);
    free(frame_->x);
    free(frame_);
    if (mtop_) free(mtop_);
    if (xcm_ != nullptr) free(xcm_);
  }

  void initAnalysis()
  {
    n_x_ = 0;

    // build molecule blocks
    if (mtop_)
    {
      // TPR path: use GROMACS molecule topology
      for (const gmx_molblock_t &molb : mtop_->molblock)
      {
        int natm_per_mol = mtop_->moltype[molb.type].atoms.nr;
        for (int i = 0; i < molb.nmol; i++) mols_.appendBlock(natm_per_mol);
        num_mol_unique_.push_back(molb.nmol);
      }
    }
    else
    {
      // Structure-file path: split by chain (each chain = one molecule instance).
      // Consecutive chains of equal atom count share a molecule type.
      if (nontpr_chain_sizes_.empty())
      {
        // No chain info (e.g. plugin did not fill chain field): one molecule.
        mols_.appendBlock(static_cast<int>(atom_names_.size()));
        num_mol_unique_.push_back(1);
      }
      else
      {
        std::size_t i = 0;
        while (i < nontpr_chain_sizes_.size())
        {
          int sz = nontpr_chain_sizes_[i];
          std::size_t j = i;
          while (j < nontpr_chain_sizes_.size() && nontpr_chain_sizes_[j] == sz) ++j;
          int nmol = static_cast<int>(j - i);
          for (int k = 0; k < nmol; k++) mols_.appendBlock(sz);
          num_mol_unique_.push_back(nmol);
          i = j;
        }
      }
    }
    nindex_ = mols_.numBlocks();

    // parse mode string (e.g. "intra+same+cross")
    printf("\nEvaluating mode selection:\n");
    std::string token;
    std::stringstream ss{ mode_ };
    while (std::getline(ss, token, '+'))
    {
      if      (token == "intra") intra_ = true;
      else if (token == "same")  same_  = true;
      else if (token == "cross") cross_ = true;
      else
      {
        printf("Wrong mode: %s\nMode must be one from: intra, same, cross. Use + to concatenate more than one, i.e. intra+cross\n", token.c_str());
        exit(1);
      }
      printf(" - found %s\n", token.c_str());
    }

    // build per-molecule lookup tables
    int mol_id = 0, molb_index = 0;
    for (int n : num_mol_unique_)
    {
      natmol2_.push_back(mols_.block(molb_index).end() - mols_.block(molb_index).begin());
      inv_num_mol_unique_.push_back(1.f / static_cast<float>(n));
      for (int j = 0; j < n; j++)
      {
        mol_id_.push_back(mol_id);
        inv_num_mol_.push_back(1.f / static_cast<float>(n));
      }
      mol_id++;
      molb_index += n;
    }

    // precompute per-type atom activity mask (skip non-backbone hydrogens)
    {
      int molb = 0, mol_first = 0;
      atom_active_.resize(natmol2_.size());
      for (std::size_t mt = 0; mt < natmol2_.size(); mt++)
      {
        atom_active_[mt].resize(natmol2_[mt]);
        const char *atomname = nullptr;
        for (int a = 0; a < natmol2_[mt]; a++)
        {
          int global_atom = mols_.block(mol_first).begin() + a;
          std::string aname;
          if (mtop_)
          {
            mtopGetAtomAndResidueName(*mtop_, global_atom, &molb, &atomname, nullptr, nullptr, nullptr);
            aname = std::string(atomname);
          }
          else
          {
            aname = atom_names_[global_atom];
          }
          const bool is_H    = (!aname.empty() && aname[0] == 'H');
          const bool is_bkbn = (aname == "H" || aname == "HN" || (!bkbn_H_.empty() && aname == bkbn_H_));
          atom_active_[mt][a] = !(is_H && !is_bkbn);
        }
        mol_first += num_mol_unique_[mt];
      }
    }

    printf("Number of different molecules %lu\n", natmol2_.size());
    bool check_same = false;
    for (std::size_t i = 0; i < natmol2_.size(); i++)
    {
      printf("mol %zu num %d size %d\n", i, num_mol_unique_[i], natmol2_[i]);
      if (num_mol_unique_[i] > 1) check_same = true;
    }
    if (!check_same && same_) same_ = false;
    if (natmol2_.size() < 2)  cross_ = false;
    if (nindex_ > 1 && (same_ || cross_))
    {
      printf("\n\n::::::::::::WARNING::::::::::::\nMore than 1 molecule found in the system.\nFix pbc before running cmdata using pbc mol\n");
      printf(":::::::::::::::::::::::::::::::\n\n");
    }
    if (same_)  std::cout << ":: activating intermat same calculations"  << std::endl;
    if (cross_) std::cout << ":: activating intermat cross calculations" << std::endl;
    if (intra_) std::cout << ":: activating intramat calculations"      << std::endl;

    // set up density bins
    n_bins_ = cmdata::indexing::n_bins(cutoff_);
    dx_     = cutoff_ / static_cast<float>(n_bins_);
    cut_sig_2_ = (cutoff_ + 0.02f) * (cutoff_ + 0.02f);
    density_bins_.resize(n_bins_);
    for (std::size_t i = 0; i < n_bins_; i++)
      density_bins_[i] = (cutoff_ / static_cast<float>(n_bins_)) * (static_cast<float>(i) + 0.5f);

    // allocate density matrices and cross index
    const std::size_t n_mol_types = natmol2_.size();
    const std::size_t n_cross     = n_mol_types * (n_mol_types - 1) / 2;
    if (same_)  { interm_same_mat_density_.resize(n_mol_types); interm_same_maxcdf_mol_.resize(n_mol_types); }
    if (cross_) { interm_cross_mat_density_.resize(n_cross);    interm_cross_maxcdf_mol_.resize(n_cross); }
    if (intra_)   intram_mat_density_.resize(n_mol_types);

    if (cross_) cross_index_.resize(n_mol_types, std::vector<int>(n_mol_types, 0));
    int cross_count = 0;
    for (std::size_t i = 0; i < n_mol_types; i++)
    {
      const auto mat_row = std::vector<std::vector<float>>(natmol2_[i], std::vector<float>(n_bins_, 0.f));
      if (same_)
      {
        interm_same_mat_density_[i].assign(natmol2_[i], mat_row);
        interm_same_maxcdf_mol_[i].assign(natmol2_[i], mat_row);
      }
      if (intra_) intram_mat_density_[i].assign(natmol2_[i], mat_row);
      for (std::size_t j = i + 1; j < n_mol_types && cross_; j++)
      {
        const auto cross_row = std::vector<std::vector<float>>(natmol2_[j], std::vector<float>(n_bins_, 0.f));
        interm_cross_mat_density_[cross_count].assign(natmol2_[i], cross_row);
        interm_cross_maxcdf_mol_[cross_count].assign(natmol2_[i], cross_row);
        cross_index_[i][j] = cross_count++;
      }
    }

    // allocate per-frame minimum-distance work buffers
    xcm_ = (rvec*)malloc(nindex_ * sizeof(rvec));
    if (same_)
    {
      frame_same_mat_.resize(n_mol_types);
      for (std::size_t i = 0; i < n_mol_types; i++)
        frame_same_mat_[i].resize(static_cast<std::size_t>(natmol2_[i]) * (natmol2_[i] + 1) / 2 * num_mol_unique_[i], 0.f);
    }
    if (cross_)
    {
      frame_cross_mat_.resize(n_cross);
      for (std::size_t i = 0; i < n_mol_types; i++)
        for (std::size_t j = i + 1; j < n_mol_types; j++)
          frame_cross_mat_[cross_index_[i][j]].resize(natmol2_[i] * natmol2_[j] * num_mol_unique_[i] * num_mol_unique_[j], 0.f);
    }

    if (!weights_path_.empty())
    {
      printf("Weights file provided. Reading weights from %s\n", weights_path_.c_str());
      weights_ = cmdata::io::read_weights_file(weights_path_);
      printf("Found %zu frame weights in file\n", weights_.size());
      float w_sum = std::accumulate(weights_.begin(), weights_.end(), 0.f);
      printf("Sum of weights amounts to %lf\n", w_sum);
    }

    printf("Finished preprocessing.\nStarting frame-by-frame analysis.\n");
  }

  void run()
  {
    // No-op callables for disabled modes — stateless lambdas, zero overhead.
    auto no_op = [](auto&&...) {};

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
      while (frame_->read_next_frame(no_pbc_, pbcType_, pbc_) == cmdata::traj::FRAME_OK)
      {
        printf("\r  Frame %d (t=%.1f ps)", frnr + 1, frame_->time);
        fflush(stdout);

        if ((frame_->time >= t_begin_ && (t_end_ < 0 || frame_->time <= t_end_)) &&
            (dt_ == 0 || std::fmod(frame_->time, dt_) == 0) &&
            (nskip_ == 0 || std::fmod(frnr, nskip_) == 0))
        {
          float weight = 1.f;
          if (!weights_.empty()) { weight = weights_[frnr]; weights_sum_ += weight; }

          // reset per-frame minimum-distance matrices
          for (auto &v : frame_same_mat_)  std::fill(v.begin(), v.end(), LARGE_DIST);
          for (auto &v : frame_cross_mat_) std::fill(v.begin(), v.end(), LARGE_DIST);

          // compute molecule centres of mass
          for (int i = 0; i < nindex_; i++)
          {
            clear_rvec(xcm_[i]);
            float tm = 0.f;
            for (int ii = mols_.block(i).begin(); ii < mols_.block(i).end(); ii++)
            {
              for (int m = 0; m < DIM; m++) xcm_[i][m] += frame_->x[ii][m];
              tm += 1.f;
            }
            for (int m = 0; m < DIM; m++) xcm_[i][m] /= tm;
          }

          // per-molecule pairwise distance computation
          #pragma omp parallel for schedule(dynamic)
          for (int i = 0; i < nindex_; i++)
            molecule_routine(i, pbc_, frame_->x, f_intra, f_same, f_cross, weight);

          // accumulate per-frame minimum distances into KDE histograms
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
      printf("\r  Read %d frames (t=%.1f ps).               \n", frnr, frame_->time);
    };

    if      ( intra_ &&  same_ &&  cross_) do_run(fi,     fs,     fc    );
    else if ( intra_ &&  same_ && !cross_) do_run(fi,     fs,     no_op );
    else if ( intra_ && !same_ &&  cross_) do_run(fi,     no_op,  fc    );
    else if ( intra_ && !same_ && !cross_) do_run(fi,     no_op,  no_op );
    else if (!intra_ &&  same_ &&  cross_) do_run(no_op,  fs,     fc    );
    else if (!intra_ &&  same_ && !cross_) do_run(no_op,  fs,     no_op );
    else if (!intra_ && !same_ &&  cross_) do_run(no_op,  no_op,  fc    );
    else                                   do_run(no_op,  no_op,  no_op );
  }

  void process_data()
  {
    std::cout << "Finished frame-by-frame analysis\n";
    std::cout << "Analyzed " << n_x_ << " frames\n";
    std::cout << "Normalizing data... " << std::endl;

    float norm = weights_.empty() ? 1.f / n_x_ : 1.f / weights_sum_;

    for (std::size_t i = 0; i < natmol2_.size(); i++)
    {
      for (int ii = 0; ii < natmol2_[i]; ii++)
      {
        for (int jj = ii; jj < natmol2_[i]; jj++)
        {
          if (same_)
          {
            cmdata::density::normalize_histo(i, ii, jj, norm, inv_num_mol_unique_[i], interm_same_maxcdf_mol_);
            cmdata::density::normalize_histo(i, ii, jj, norm, 1.f, interm_same_mat_density_);
            float sum = 0.f;
            for (std::size_t k = 0; k < n_bins_; k++)
            {
              sum += dx_ * interm_same_maxcdf_mol_[i][ii][jj][k];
              if (sum > 1.f) sum = 1.f;
              interm_same_maxcdf_mol_[i][ii][jj][k] = sum;
            }
            interm_same_mat_density_[i][jj][ii] = interm_same_mat_density_[i][ii][jj];
            interm_same_maxcdf_mol_[i][jj][ii]  = interm_same_maxcdf_mol_[i][ii][jj];
          }
          if (intra_)
          {
            cmdata::density::normalize_histo(i, ii, jj, norm, 1.f, intram_mat_density_);
            intram_mat_density_[i][jj][ii] = intram_mat_density_[i][ii][jj];
          }
        }
      }
      if (cross_)
      {
        for (std::size_t j = i + 1; j < natmol2_.size(); j++)
        {
          for (int ii = 0; ii < natmol2_[i]; ii++)
          {
            for (int jj = 0; jj < natmol2_[j]; jj++)
            {
              cmdata::density::normalize_histo(cross_index_[i][j], ii, jj, norm, 1.f, interm_cross_mat_density_);
              cmdata::density::normalize_histo(cross_index_[i][j], ii, jj, norm, inv_num_mol_unique_[i], interm_cross_maxcdf_mol_);
              float sum = 0.f;
              for (std::size_t k = 0; k < n_bins_; k++)
              {
                sum += dx_ * interm_cross_maxcdf_mol_[cross_index_[i][j]][ii][jj][k];
                if (sum > 1.f) sum = 1.f;
                interm_cross_maxcdf_mol_[cross_index_[i][j]][ii][jj][k] = sum;
              }
            }
          }
        }
      }
    }
  }

  void write_output(const std::string &output_prefix)
  {
    std::cout << "Writing data... " << std::endl;

    // select write functions (text or HDF5) once up front
    auto* wf_intra = &cmdata::io::f_write_intra;
    auto* wf_same  = &cmdata::io::f_write_inter_same;
    auto* wf_cross = &cmdata::io::f_write_inter_cross;
    #ifdef USE_HDF5
    if (h5_)
    {
      wf_intra = &cmdata::io::f_write_intra_HDF5;
      wf_same  = &cmdata::io::f_write_inter_same_HDF5;
      wf_cross = &cmdata::io::f_write_inter_cross_HDF5;
    }
    #endif

    for (std::size_t i = 0; i < natmol2_.size(); i++)
    {
      std::cout << "Molecule " << i << ":" << std::endl;
      float progress = 0.f;
      cmdata::io::print_progress_bar(0.f);
      for (int ii = 0; ii < natmol2_[i]; ii++)
      {
        float new_progress = static_cast<float>(ii + 1) / static_cast<float>(natmol2_[i]);
        if (new_progress - progress > 0.01f) { progress = new_progress; cmdata::io::print_progress_bar(progress); }
        if (intra_) wf_intra(output_prefix, i, ii, density_bins_, natmol2_, intram_mat_density_);
        if (same_)  wf_same (output_prefix, i, ii, density_bins_, natmol2_, interm_same_mat_density_, interm_same_maxcdf_mol_);
      }
      if (cross_)
      {
        for (std::size_t j = i + 1; j < natmol2_.size(); j++)
          for (int ii = 0; ii < natmol2_[i]; ii++)
            wf_cross(output_prefix, i, j, ii, density_bins_, natmol2_, cross_index_, interm_cross_mat_density_, interm_cross_maxcdf_mol_);
      }
      cmdata::io::print_progress_bar(1.f);  // closes the bar line with \n
    }
    std::cout << "Finished!" << std::endl;
  }
};

} // namespace cmdata

#endif // _CMDATA_CMDATA_HPP
