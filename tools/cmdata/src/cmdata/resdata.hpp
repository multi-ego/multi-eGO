#ifndef _RESDATA_CMDATA_HPP
#define _RESDATA_CMDATA_HPP

// gromacs includes
#include <gromacs/trajectoryanalysis/topologyinformation.h>
#include <gromacs/math/vec.h>
#include <gromacs/pbcutil/pbc.h>
#include <gromacs/fileio/tpxio.h>

// cmdata includes
#include "io.hpp"
#include "indexing.hpp"
#include "parallel.hpp"
#include "density.hpp"
#include "mindist.hpp"
#include "xtc_frame.hpp"
#include "function_types.hpp"

// standard library imports
#include <iostream>
#include <omp.h>
#include <thread>
#include <mutex>
#include <condition_variable>
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

void prob_same(
  int i, std::size_t a_i, std::size_t a_j, double dx2, double weight,
  const std::vector<int> &mol_id_, const std::vector<int> &nresmol2_,
  std::vector<std::vector<std::mutex>> &frame_same_mutex_,
  std::vector<std::vector<std::vector<double>>> &interm_prob_same_,
  const std::vector<std::vector<std::vector<int>>> &interm_frame_same_,
  const std::vector<std::vector<int>> &rndx_
)
{

  std::size_t same_mutex_index = cmdata::indexing::mutex_access(mol_id_[i], r_i, r_j, nresmol2_);
  std::unique_lock lock(frame_same_mutex_[mol_id_[i]][same_mutex_index]);
  interm_prob_same_[mol_id_[i]][r_i][r_j]++;
  lock.unlock();
}

void prob_cross(
  int i, int j, std::size_t a_i, std::size_t a_j, double dx2, double weight,
  const std::vector<int> &mol_id_, const std::vector<int> &nresmol2_,
  const std::vector<std::vector<int>> &cross_index_,
  std::vector<std::vector<std::mutex>> &frame_cross_mutex_,
  std::vector<std::vector<std::vector<double>>> &interm_cross_mindist_,
  const std::vector<std::vector<int>> &rndx_
)
{
  std::size_t r_i = rndx_[mol_id_[i]][a_i];
  std::size_t r_j = rndx_[mol_id_[j]][a_j];
  std::size_t cross_mutex_index = cmdata::indexing::mutex_access(mol_id_[j], a_i, a_j, nresmol2_);
  std::unique_lock lock(frame_cross_mutex_[cross_index_[mol_id_[i]][mol_id_[j]]][cross_mutex_index]);
  interm_cross_mindist_[cross_index_[mol_id_[i]][mol_id_[j]]][r_i][r_j]++;
  lock.unlock();
}

class ResData
{
private:
  // general fields
  int n_x_;
  float dt_, t_begin_, t_end_;
  int nskip_;
  gmx_mtop_t *mtop_;
  rvec *xcm_ = nullptr;
  double cutoff_, mol_cutoff_, mcut2_, cut_sig_2_;

  // molecule number fields
  int nindex_;
  // std::vector<uint> selection_;
  gmx::RangePartitioning mols_;
  std::vector<int> natmol2_;
  std::vector<int> nresmol2_;
  std::vector<int> mol_id_;
  std::vector<int> num_mol_unique_;
  std::vector<double> inv_num_mol_unique_;
  std::vector<double> inv_num_mol_;
  
  std::vector<std::vector<int>> rndx_;
  std::vector<std::vector<int>> res_sizes_;

  // weights fields
  double weights_sum_;
  std::string weights_path_;
  std::vector<double> weights_;

  // pbc fields
  bool no_pbc_;
  t_pbc *pbc_;
  PbcType pbcType_;

  // frame fields
  cmdata::xtc::Frame *frame_;
  XDRFILE *trj_;

  // density fields
  double dx_;
  std::size_t n_bins_;
  std::vector<double> density_bins_;
  std::vector<std::vector<int>> cross_index_;

  using cmdata_matrix = std::vector<std::vector<std::vector<double>>>;
  cmdata_matrix interm_same_mat_density_;
  cmdata_matrix interm_cross_mat_density_;
  // cmdata_matrix intram_mat_density_;
  // cmdata_matrix interm_same_maxcdf_mol_;
  // cmdata_matrix interm_cross_maxcdf_mol_;

  // temporary containers for maxcdf operations
  std::vector<std::vector<std::vector<int>>> frame_same_mat_;
  std::vector<std::vector<std::vector<int>>> frame_cross_mat_;
  std::vector<std::vector<std::mutex>> frame_same_mutex_;
  std::vector<std::vector<std::mutex>> frame_cross_mutex_;

  // parallelization fields
  int num_threads_;
  int num_mol_threads_;
  std::vector<std::thread> threads_;
  std::vector<std::thread> mol_threads_;
  cmdata::parallel::Semaphore semaphore_;
  // std::vector<cmdata::indexing::SameThreadIndices> same_thread_indices_;
  // std::vector<cmdata::indexing::CrossThreadIndices> cross_thread_indices_;

  // mode selection, booleans and functions
  std::string mode_;
  bool intra_ = false, same_ = false, cross_ = false;

  // function types
  using ftype_intra_ = cmdata::ftypes::function_traits<decltype(&cmdata::density::intra_mol_routine)>;
  using ftype_same_ = cmdata::ftypes::function_traits<decltype(&prob_same)>;
  using ftype_cross_ = cmdata::ftypes::function_traits<decltype(&prob_cross)>;

  std::function<ftype_intra_::signature> f_intra_mol_;
  std::function<ftype_same_::signature> f_inter_mol_same_;
  std::function<ftype_cross_::signature> f_inter_mol_cross_;

  static void molecule_routine(
    const int i, const int nindex_, t_pbc *pbc, rvec *x, const std::vector<double> &inv_num_mol_, const double cut_sig_2_, 
    const std::vector<int> &natmol2_, const std::vector<int> &num_mol_unique_, const std::vector<int> &mol_id_, 
    const std::vector<std::vector<int>> &cross_index_, const std::vector<double> &density_bins_, const double mcut2_, 
    rvec *xcm_, const gmx::RangePartitioning &mols_, gmx_mtop_t *mtop_, 
    // std::vector<std::vector<double>> &frame_same_mat_, 
    std::vector<std::vector<std::mutex>> &frame_same_mutex_,
    // cmdata_matrix &intram_mat_density_, 
    cmdata_matrix &interm_same_mat_density_, 
    // std::vector<std::vector<double>> &frame_cross_mat_,
    std::vector<std::vector<std::mutex>> &frame_cross_mutex_, cmdata_matrix &interm_cross_mat_density_, cmdata::parallel::Semaphore &semaphore_,
    const std::function<ftype_intra_::signature> &f_intra_mol_, const std::function<ftype_same_::signature> &f_inter_mol_same_,
    const std::function<ftype_cross_::signature> &f_inter_mol_cross_, double weight, 
    const std::vector<std::vector<int>> &rndx_, const std::vector<std::vector<int>> &res_sizes_
  )
  {
    semaphore_.acquire();
    const char * atomname;
    int tmp_i = 0;
    std::size_t mol_i = i, mol_j = 0;
    while ( static_cast<int>(mol_i) - num_mol_unique_[tmp_i] >= 0  )
    {
      mol_i -= num_mol_unique_[tmp_i];
      tmp_i++;
      if (tmp_i == num_mol_unique_.size()) break;
    }
    if (mol_i == num_mol_unique_[mol_id_[i]]) mol_i = 0;
    int molb = 0;
    int resi_id = 0;
    int resj_id = 0;
    int resi_inc = res_sizes_[mol_id_[i]][0];
    /* Loop over molecules  */
    for (int j = 0; j < nindex_; j++)
    {
      int resi_inc = res_sizes_[mol_id_[j]][0];
      if (j!=0)
        if (mol_j == num_mol_unique_[mol_id_[j-1]]) mol_j = 0;

      /* intermolecular interactions are evaluated only among neighbour molecules */
      if (i!=j)
      {
        rvec dx;
        if (pbc != nullptr) pbc_dx(pbc, xcm_[i], xcm_[j], dx);
        else rvec_sub(xcm_[i], xcm_[j], dx);
        double dx2 = iprod(dx, dx);
        if (dx2 > mcut2_) continue;
      }
      /* for molecules of different specie we fill half a matrix */
      if (mol_id_[i] != mol_id_[j] && j < i) continue;
      // std::size_t a_i = 0;
      // GMX_RELEASE_ASSERT(mols_.numBlocks() > 0, "Cannot access index[] from empty mols");

      /* cycle over the atoms of a molecule i */
      for (std::size_t ii = mols_.block(i).begin(); ii < mols_.block(i).end(); ii++)
      {
        // std::size_t a_j = 0;
        // mtopGetAtomAndResidueName(*mtop_, ii, &molb, &atomname, nullptr, nullptr, nullptr);
        // if (atomname[0] == 'H')
        // {
          // a_i++;
          // continue;
        // }
        /* cycle over the atoms of a molecule j */
        for (std::size_t jj = mols_.block(j).begin(); jj < mols_.block(j).end(); jj++)
        {
          // mtopGetAtomAndResidueName(*mtop_, jj, &molb, &atomname, nullptr, nullptr, nullptr);
          // if (atomname[0] == 'H')
          // {
            // a_j++;
            // continue;
          // }
          // std::size_t delta  = a_i - a_j;
          rvec sym_dx;
          if (pbc != nullptr) pbc_dx(pbc, x[ii], x[jj], sym_dx);
          else rvec_sub(x[ii], x[jj], sym_dx);
          double dx2 = iprod(sym_dx, sym_dx);
          if (i==j) 
          {
            if (dx2 < cut_sig_2_)
            { // intra molecule species
              // f_intra_mol_(i, a_i, a_j, dx2, weight, mol_id_, natmol2_, density_bins_, inv_num_mol_, frame_same_mutex_, intram_mat_density_);
            }
          }
          else
          {
            if (mol_id_[i]==mol_id_[j])
            { // inter same molecule specie
              if (dx2 < cut_sig_2_)
              {
                f_inter_mol_same_(i, a_i, a_j, dx2, weight, mol_id_, natmol2_, frame_same_mutex_, interm_same_mat_density_, frame_same_mat_, rndx_);
              }
              // if (delta!=0.) {
              //   // this is to account for inversion atom/molecule
              //   if (pbc != nullptr) pbc_dx(pbc, x[ii-delta], x[jj+delta], sym_dx);
              //   else rvec_sub(x[ii-delta], x[jj+delta], sym_dx);
              //   dx2 = iprod(sym_dx, sym_dx);
              //   if (dx2 < cut_sig_2_)
              //   {
              //     // f_inter_mol_same_(
              //     //   i, mol_i, a_i, a_j, dx2, weight, mol_id_, natmol2_, density_bins_, frame_same_mutex_, frame_same_mat_, interm_same_mat_density_
              //     // );
              //     f_inter_mol_same_(i, a_i, a_j, dx2, weight, mol_id_, natmol2_, frame_same_mutex_, interm_same_mat_density_, frame_same_mat_, rndx_);
              //   }
              // }
            } 
            else
            { // inter cross molecule species
              if (dx2 < cut_sig_2_)
              {
                // f_inter_mol_cross_(i, j, a_i, a_j, dx2, weight, mol_id_, natmol2_, cross_index_, frame_cross_mutex_, interm_cross_mat_density_, rndx_);
              }
            }
          }
          // ++a_j;
        }
        // ++a_i;
      }
      ++mol_j;
    }
    semaphore_.release();
  }

  // stolen from carlo's code
  std::vector<int> res_ndx(t_atoms& atoms)
  {
    std::vector<int> rndx;
    int  i, r0;

    if (atoms.nr <= 0)
    {
        return std::vector<int>();
    }
    r0 = atoms.atom[0].resind;
    for (i = 0; (i < atoms.nr); i++)
    {
      rndx.push_back(atoms.atom[i].resind - r0);
    }

    return rndx;
  }

static std::vector<int> res_natm(t_atoms &atoms)
{
  // int* natm;
  std::vector<int> natm(atoms.nres, 0);
  int  i, j, r0;
  
  if (atoms.nr <= 0)
  {
    return natm;
  }
  // snew(natm, atoms->nres;);
  r0 = atoms.atom[0].resind;
  j  = 0;
  for (i = 0; (i < atoms.nres); i++)
  {
    while ((atoms.atom[j].resind) - r0 == i)
    {
      natm[i]++;
      j++;
    }
  }

  return natm;
}

public:
  ResData(
    const std::string &top_path, const std::string &traj_path,
    double cutoff, double mol_cutoff, int nskip, int num_threads, int num_mol_threads,
    int dt, const std::string &mode, const std::string &weights_path, 
    bool no_pbc, float t_begin, float t_end
  ) : cutoff_(cutoff), mol_cutoff_(mol_cutoff), nskip_(nskip), num_threads_(num_threads), num_mol_threads_(num_mol_threads),
      mode_(mode), weights_path_(weights_path),
      no_pbc_(no_pbc), dt_(dt), t_begin_(t_begin), t_end_(t_end)
  {
    bool bTop_;
    matrix boxtop_;
    mtop_ = (gmx_mtop_t*)malloc(sizeof(gmx_mtop_t));
    TpxFileHeader header = readTpxHeader(top_path.c_str(), true);
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
    frame_->nframe = nframe;
    frame_->offsets = offsets;

    trj_ = xdrfile_open(traj_path.c_str(), "r");
    initAnalysis();
  }

  ~ResData()
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
      rndx_.push_back(res_ndx(mtop_->moltype[molb.type].atoms));
      res_sizes_.push_back(res_natm(mtop_->moltype[molb.type].atoms));
    }

    // number of molecules
    nindex_ = mols_.numBlocks();

    if (num_threads_ > std::thread::hardware_concurrency())
    {
      num_threads_ = std::thread::hardware_concurrency();
      std::cout << "Maximum thread number surpassed. Scaling num_threads down to " << num_threads_ << std::endl;
    }
    if (num_mol_threads_ > std::thread::hardware_concurrency())
    {
      num_mol_threads_ = std::thread::hardware_concurrency();
      std::cout << "Maximum thread number surpassed. Scaling num_mol_threads down to " << num_mol_threads_ << std::endl;
    }
    if (num_mol_threads_ > nindex_)
    {
      num_mol_threads_ = nindex_;
      std::cout << "Number of molecule threads surpassed number of molecules. Setting num_mol_threads to " << num_mol_threads_ << std::endl;
    }
    threads_.resize(num_threads_);
    mol_threads_.resize(nindex_);
    semaphore_.set_counter(num_mol_threads_);
    std::cout << "Using " << num_threads_ << " threads and " << num_mol_threads_ << " molecule threads" << std::endl;

    // set up mode selection
    f_intra_mol_ = cmdata::ftypes::do_nothing<ftype_intra_>();
    f_inter_mol_same_ = cmdata::ftypes::do_nothing<ftype_same_>();
    f_inter_mol_cross_ = cmdata::ftypes::do_nothing<ftype_cross_>();

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
      nresmol2_.push_back(rndx_[mol_id].back() + 1);
      inv_num_mol_unique_.push_back(1. / static_cast<double>(i));
      for ( int j = 0; j < i; j++ )
      {
        mol_id_.push_back(mol_id);
        inv_num_mol_.push_back(1. / static_cast<double>(i));
      }
      mol_id++;
      molb_index += i;
    }

    printf("Number of different molecules %lu\n", natmol2_.size());
    bool check_same = false;
    for(std::size_t i=0; i<natmol2_.size();i++) {
      printf("mol %lu num %u size %u\n", i, num_mol_unique_[i], natmol2_[i]);
      if(num_mol_unique_[i]>1) check_same = true;
    }
    if(!check_same && same_) same_ = false;
    if(natmol2_.size()<2) cross_ = false; 

    if (same_)
    {
      f_inter_mol_same_ = prob_same;
      std::cout << " :: activating intermat same calculations" << std::endl;
      interm_same_mat_density_.resize(natmol2_.size());
      // interm_same_maxcdf_mol_.resize(natmol2_.size());
    }
    if (cross_)
    {
      f_inter_mol_cross_ = prob_cross;
      std::cout << " :: activating intermat cross calculations" << std::endl;
      interm_cross_mat_density_.resize((natmol2_.size() * (natmol2_.size() - 1)) / 2);
      // interm_cross_maxcdf_mol_.resize((natmol2_.size() * (natmol2_.size() - 1)) / 2);
    }
    if (intra_) 
    {
      // f_intra_mol_ = cmdata::density::intra_mol_routine;
      // std::cout << " :: activating intramat calculations" << std::endl;
      // intram_mat_density_.resize(natmol2_.size());
    }

    density_bins_.resize(cmdata::indexing::n_bins(cutoff_));
    for (std::size_t i = 0; i < density_bins_.size(); i++)
      density_bins_[i] = cutoff_ / static_cast<double>(density_bins_.size()) * static_cast<double>(i) + cutoff_ / static_cast<double>(density_bins_.size() * 2);

    int cross_count = 0;
    if (cross_) cross_index_.resize(natmol2_.size(), std::vector<int>(natmol2_.size(), 0));
    for ( std::size_t i = 0; i < natmol2_.size(); i++ )
    {
      if (same_)
      {
        // std::vector<std::vector<std::vector<int>>>
        interm_same_mat_density_[i].resize(nresmol2_[i], std::vector<double>(nresmol2_[i], 0));
        // interm_same_maxcdf_mol_[i].resize(natmol2_[i], std::vector<std::vector<double>>(natmol2_[i], std::vector<double>(cmdata::indexing::n_bins(cutoff_), 0)));
      }
      // if (intra_) intram_mat_density_[i].resize(natmol2_[i], std::vector<std::vector<double>>(natmol2_[i], std::vector<double>(cmdata::indexing::n_bins(cutoff_), 0)));
      for ( std::size_t j = i + 1; j < natmol2_.size() && cross_; j++ )
      {
        interm_cross_mat_density_[cross_count].resize(nresmol2_[i], std::vector<double>(nresmol2_[j], 0));
        // interm_cross_mat_density_[cross_count].resize(natmol2_[i], std::vector<std::vector<double>>(natmol2_[j], std::vector<double>(cmdata::indexing::n_bins(cutoff_), 0)));
        // interm_cross_maxcdf_mol_[cross_count].resize(natmol2_[i], std::vector<std::vector<double>>(natmol2_[j], std::vector<double>(cmdata::indexing::n_bins(cutoff_), 0)));
        cross_index_[i][j] = cross_count;
        cross_count++;
      }
    }

    n_bins_ = cmdata::indexing::n_bins(cutoff_);
    dx_ = cutoff_ / static_cast<double>(n_bins_);

    mcut2_ = mol_cutoff_ * mol_cutoff_;
    cut_sig_2_ = (cutoff_ + 0.02) * (cutoff_ + 0.02);
    xcm_ = (rvec*)malloc(nindex_ * sizeof(rvec));

    if (same_) frame_same_mat_.resize(natmol2_.size());
    if (intra_ || same_) frame_same_mutex_.resize(nresmol2_.size());
    // if (cross_) frame_cross_mat_.resize(cross_index_.size());
    if (cross_) frame_cross_mutex_.resize(cross_index_.size());
    std::size_t sum_cross_mol_sizes = 0;
    for ( std::size_t i = 0; i < natmol2_.size(); i++ )
    {
      // if (same_) frame_same_mat_[i].resize(natmol2_[i] *  natmol2_[i] * num_mol_unique_[i], 0);
      if (intra_ || same_) frame_same_mutex_[i] = std::vector<std::mutex>(nresmol2_[i] *  nresmol2_[i]);
      for ( std::size_t j = i+1; j < nresmol2_.size() && cross_; j++ )
      {
        // frame_cross_mat_[cross_index_[i][j]].resize(natmol2_[i] * natmol2_[j] * num_mol_unique_[i] * num_mol_unique_[j], 0);
        frame_cross_mutex_[cross_index_[i][j]] = std::vector<std::mutex>(nresmol2_[i] * nresmol2_[j]);
      }
    }

    if (weights_path_ != "")
    {
      printf("Weights file provided. Reading weights from %s\n", weights_path_.c_str());
      weights_ = cmdata::io::read_weights_file(weights_path_);
      printf("Found %li frame weights in file\n", weights_.size());
      double w_sum = std::accumulate(std::begin(weights_), std::end(weights_), 0.0, std::plus<>());
      printf("Sum of weights amounts to %lf\n", w_sum);
      weights_sum_ = 0.;
    }

    std::cout << "Calculating threading indices" << std::endl;
    /* calculate the mindist accumulation indices */
    std::size_t num_ops_same = 0;
    for ( std::size_t im = 0; im < natmol2_.size(); im++ ) num_ops_same += num_mol_unique_[im] * ( nresmol2_[im] * ( nresmol2_[im] + 1 ) ) / 2;
    int n_per_thread_same = (same_) ? num_ops_same / num_threads_  : 0;
    int n_threads_same_uneven = (same_) ? num_ops_same % num_threads_ : 0;
    std::size_t start_mti_same = 0, start_im_same = 0, end_mti_same = 1, end_im_same = 1; 
    std::size_t start_i_same = 0, start_j_same = 0, end_i_same = 0, end_j_same = 0;
    int num_ops_cross = 0;
    for ( std::size_t im = 0; im < natmol2_.size(); im++ )
    {
      for ( std::size_t jm = im + 1; jm < natmol2_.size(); jm++ )
      {
        num_ops_cross += num_mol_unique_[im] * nresmol2_[im] * num_mol_unique_[jm] * nresmol2_[jm];
      }
    }
    int n_per_thread_cross = (cross_) ? num_ops_cross / num_threads_ : 0;
    int n_threads_cross_uneven = (cross_) ? num_ops_cross % num_threads_ : 0;

    std::size_t start_mti_cross = 0, start_mtj_cross = 1, start_im_cross = 0, start_jm_cross = 0, start_i_cross = 0, start_j_cross = 0;
    std::size_t end_mti_cross = 1, end_mtj_cross = 2, end_im_cross = 1, end_jm_cross = 1, end_i_cross = 0, end_j_cross = 0;

    // for ( int tid = 0; tid < num_threads_; tid++ )
    // {
    //   /* calculate same indices */
    //   int n_loop_operations_same = n_per_thread_same + (tid < n_threads_same_uneven ? 1 : 0);
    //   long int n_loop_operations_total_same = n_loop_operations_same;
    //   while ( natmol2_[end_mti_same - 1] - static_cast<int>(end_j_same) <= n_loop_operations_same )
    //   {
    //     int sub_same = natmol2_[end_mti_same - 1] - static_cast<int>(end_j_same);
    //     n_loop_operations_same -= sub_same;
    //     end_i_same++;
    //     end_j_same = end_i_same;
    //     if (static_cast<int>(end_j_same) == natmol2_[end_mti_same - 1])
    //     {
    //       end_im_same++;
    //       end_i_same = 0;
    //       end_j_same = 0;
    //     }
    //     if (static_cast<int>(end_im_same) -1 == num_mol_unique_[end_mti_same - 1])
    //     {
    //       end_mti_same++;
    //       end_im_same = 1;
    //       end_i_same = 0;
    //       end_j_same = 0;
    //     }
    //     if (n_loop_operations_same == 0) break;
    //   }
    //   end_j_same += n_loop_operations_same;  
    //   /* calculate cross indices */
    //   int n_loop_operations_total_cross = n_per_thread_cross + ( tid < n_threads_cross_uneven ? 1 : 0 );
    //   if (natmol2_.size() > 1)
    //   {
    //     int n_loop_operations_cross = n_loop_operations_total_cross;
    //     while ( natmol2_[end_mti_cross-1] * natmol2_[end_mtj_cross-1] - (natmol2_[end_mtj_cross-1] * static_cast<int>(end_i_cross) + static_cast<int>(end_j_cross)) <= n_loop_operations_cross )
    //     {
    //       int sub_cross = natmol2_[end_mti_cross-1] * natmol2_[end_mtj_cross-1] - (natmol2_[end_mtj_cross-1] * static_cast<int>(end_i_cross) + static_cast<int>(end_j_cross));
    //       n_loop_operations_cross -= sub_cross;

    //       end_jm_cross++;
    //       end_i_cross = 0;
    //       end_j_cross = 0;

    //       // case jm is above max
    //       if (end_jm_cross > num_mol_unique_[end_mtj_cross - 1])
    //       {
    //         // end_im_cross++;
    //         end_mtj_cross++;
    //         end_jm_cross = 1;
    //         end_i_cross = 0;
    //         end_j_cross = 0;
    //       }
    //       if (end_mtj_cross > natmol2_.size())
    //       {
    //         end_im_cross++;
    //         end_mtj_cross = end_mti_cross + 1;
    //         end_jm_cross = 1;
    //         end_i_cross = 0;
    //         end_j_cross = 0;
    //       }
    //       if (end_im_cross > num_mol_unique_[end_mti_cross - 1])
    //       {
    //         end_mti_cross++;
    //         end_mtj_cross = end_mti_cross + 1;
    //         end_im_cross = 1;
    //         end_jm_cross = 1;
    //         end_i_cross = 0;
    //         end_j_cross = 0;
    //       }
    //       if (end_mti_cross == natmol2_.size()) break;
    //       if (n_loop_operations_cross == 0) break;
    //     }

    //     // calculate overhangs and add them
    //     if (end_mti_cross < natmol2_.size())
    //     {
    //       end_i_cross += n_loop_operations_cross / natmol2_[end_mtj_cross-1];
    //       end_j_cross += n_loop_operations_cross % natmol2_[end_mtj_cross-1];
    //       end_i_cross += end_j_cross / natmol2_[end_mtj_cross-1];
    //       end_j_cross %= natmol2_[end_mtj_cross-1];
    //     }
    //   }

    //   if (same_)
    //   {
    //     same_thread_indices_.push_back({
    //       start_mti_same, start_im_same, start_i_same, start_j_same, end_mti_same,
    //       end_im_same, end_i_same, end_j_same, n_loop_operations_total_same
    //     });
    //   }
    //   else
    //   {
    //     same_thread_indices_.push_back({
    //       0, 0, 0, 0, 0, 0, 0, 0, 0
    //     });
    //   }
    //   if (cross_ && natmol2_.size() > 1)
    //   {
    //     cross_thread_indices_.push_back({
    //       start_mti_cross, start_mtj_cross, start_im_cross, start_jm_cross, start_i_cross,
    //       start_j_cross, end_mti_cross, end_mtj_cross, end_im_cross, end_jm_cross, end_i_cross,
    //       end_j_cross, n_loop_operations_total_cross
    //     });
    //   }
    //   else
    //   {
    //     cross_thread_indices_.push_back({
    //       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    //     });
    //   }

    //   /* set new starts */
    //   start_mti_same = end_mti_same - 1;
    //   start_im_same = end_im_same - 1;
    //   start_i_same = end_i_same;
    //   start_j_same = end_j_same;

    //   start_mti_cross = end_mti_cross - 1;
    //   start_mtj_cross = end_mtj_cross - 1;
    //   start_im_cross = end_im_cross - 1;
    //   start_jm_cross = end_jm_cross - 1;
    //   start_i_cross = end_i_cross;
    //   start_j_cross = end_j_cross;
    // }

    printf("Finished preprocessing. Starting frame-by-frame analysis.\n");
  }

  void run()
  {
    std::cout << "Running frame-by-frame analysis" << std::endl;
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
      if ((frame_->time >= t_begin_ && (t_end_ < 0 || frame_->time <= t_end_ )) && // within time borders 
          ( dt_ == 0 || std::fmod(frame_->time, dt_) == 0) && (nskip_ == 0 || std::fmod(frnr, nskip_) == 0)) // skip frames
      {
        double weight = 1.0;
        if (!weights_.empty())
        {
          weight = weights_[frnr];
          weights_sum_ += weight;
        }
        /* resetting the per frame vector to zero */
        // for ( std::size_t i = 0; i < frame_same_mat_.size(); i++ )
        // {
        //   #pragma omp parallel for num_threads(std::min(num_threads_, static_cast<int>(frame_same_mat_[i].size())))
        //   for ( std::size_t j = 0; j < frame_same_mat_[i].size(); j++ ) frame_same_mat_[i][j] = 100.;
        // }
        // for ( std::size_t i = 0; i < frame_cross_mat_.size(); i++ )
        // {
        //   #pragma omp parallel for num_threads(std::min(num_threads_, static_cast<int>(frame_cross_mat_[i].size())))
        //   for ( std::size_t j = 0; j < frame_cross_mat_[i].size(); j++ ) frame_cross_mat_[i][j] = 100.;
        // }
        #pragma omp parallel for num_threads(std::min(num_threads_, nindex_))
        for (int i = 0; i < nindex_; i++)
        {
          clear_rvec(xcm_[i]);
          double tm = 0.;
          for (int ii = mols_.block(i).begin(); ii < mols_.block(i).end(); ii++)
          {
            for (int m = 0; (m < DIM); m++)
            {
              xcm_[i][m] += frame_->x[ii][m];
            }
            tm += 1.0;
          }
          for (int m = 0; (m < DIM); m++)
          {
            xcm_[i][m] /= tm;
          }
        }
        /* start loop for each molecule */
        for (int i = 0; i < nindex_; i++)
        {
          /* start molecule thread*/
          // mol_threads_[i] = std::thread(molecule_routine, i, nindex_, pbc_, frame_->x, std::cref(inv_num_mol_), 
          // cut_sig_2_, std::cref(natmol2_), std::cref(num_mol_unique_), std::cref(mol_id_), std::cref(cross_index_),
          // std::cref(density_bins_), mcut2_, xcm_, mols_, mtop_, std::ref(frame_same_mat_),
          // std::ref(frame_same_mutex_), std::ref(intram_mat_density_), std::ref(interm_same_mat_density_), std::ref(frame_cross_mat_), 
          // std::ref(frame_cross_mutex_), std::ref(interm_cross_mat_density_), std::ref(semaphore_), std::cref(f_intra_mol_),
          // std::cref(f_inter_mol_same_), std::cref(f_inter_mol_cross_), weight);
    
    // molecule routine parameters
    // const int i, const int nindex_, t_pbc *pbc, rvec *x, const std::vector<double> &inv_num_mol_, const double cut_sig_2_, 
    // const std::vector<int> &natmol2_, const std::vector<int> &num_mol_unique_, const std::vector<int> &mol_id_, 
    // const std::vector<std::vector<int>> &cross_index_, const std::vector<double> &density_bins_, const double mcut2_, 
    // rvec *xcm_, const gmx::RangePartitioning &mols_, gmx_mtop_t *mtop_, 
    // std::vector<std::vector<std::mutex>> &frame_same_mutex_,
    // cmdata_matrix &interm_same_mat_density_, 
    // std::vector<std::vector<std::mutex>> &frame_cross_mutex_, cmdata_matrix &interm_cross_mat_density_, cmdata::parallel::Semaphore &semaphore_,
    // const std::function<ftype_intra_::signature> &f_intra_mol_, const std::function<ftype_same_::signature> &f_inter_mol_same_,
    // const std::function<ftype_cross_::signature> &f_inter_mol_cross_, double weight, std::vector<std::vector<int>> rndx_
          molecule_routine(
            i, nindex_, pbc_, frame_->x, std::cref(inv_num_mol_), cut_sig_2_, std::cref(nresmol2_), std::cref(num_mol_unique_), 
            std::cref(mol_id_), std::cref(cross_index_), std::cref(density_bins_), mcut2_, xcm_, mols_, mtop_, 
            std::ref(frame_same_mutex_), std::ref(interm_same_mat_density_), 
            std::ref(frame_cross_mutex_), std::ref(interm_cross_mat_density_), std::ref(semaphore_), 
            std::cref(f_intra_mol_), std::cref(f_inter_mol_same_), std::cref(f_inter_mol_cross_), weight, rndx_
          );
          /* end molecule thread */
        }
        /* join molecule threads */
        // for ( auto &thread : mol_threads_ ) thread.join();

        /* calculate the mindist accumulation indices */
        // for ( int tid = 0; tid < num_threads_; tid++ )
        // {
        //   threads_[tid] = std::thread(
        //     cmdata::mindist::mindist_kernel, std::cref(same_thread_indices_[tid]), std::cref(cross_thread_indices_[tid]),
        //     weight, std::cref(natmol2_), std::cref(density_bins_), std::cref(num_mol_unique_),
        //     std::cref(frame_same_mat_), std::ref(frame_same_mutex_), std::ref(interm_same_maxcdf_mol_),
        //     std::cref(cross_index_), std::cref(frame_cross_mat_), std::ref(frame_cross_mutex_), std::ref(interm_cross_maxcdf_mol_)
        //   );
        // }
        // for ( auto &thread : threads_ ) thread.join();
        ++n_x_;
      }
      ++frnr;
    }

    cmdata::io::print_progress_bar(1.0);
  }

  static void normalize_resdata(
    std::size_t i, int ii, double norm, double inv_num_mol_same,
    std::vector<std::vector<std::vector<double>>> &data
  )
  {
    // std::transform(
    //   std::begin(data[i][ii]), std::end(data[i][ii]), std::begin(data[i][ii]), 
    //   [norm, inv_num_mol_same]( auto &c ){ return c * norm * inv_num_mol_same; }
    // );
  }

  void process_data()
  {
    std::cout << "Finished frame-by-frame analysis\n";
    std::cout << "Analyzed " << n_x_ << " frames\n";
    std::cout << "Normalizing data... " << std::endl;
    // normalisations
    double norm = ( weights_.empty() ) ? 1. / n_x_ : 1. / weights_sum_;

    using ftype_norm = cmdata::ftypes::function_traits<decltype(&normalize_resdata)>;
    std::function<ftype_norm::signature> f_empty = cmdata::ftypes::do_nothing<ftype_norm>();

    std::function<ftype_norm::signature> normalize_intra = (intra_) ? normalize_resdata : f_empty;
    std::function<ftype_norm::signature> normalize_inter_same = (same_) ? normalize_resdata : f_empty;
    std::function<ftype_norm::signature> normalize_inter_cross = (cross_) ? normalize_resdata : f_empty;

    for (std::size_t i = 0; i < nresmol2_.size(); i++)
    {
      for (int ii = 0; ii < nresmol2_[i]; ii++)
      {
        // for (int jj = ii; jj < nresmol2_[i]; jj++)
        // {
          double inv_num_mol_same = inv_num_mol_unique_[i];
          // normalize_inter_same(i, ii, norm, inv_num_mol_same, interm_same_maxcdf_mol_);
          // normalize_inter_same(i, ii, norm, inv_num_mol_same, interm_same_mat_density_);
          // normalize_intra(i, ii, norm, 1.0, intram_mat_density_);

          // for ()



          double sum = 0.0;
          for ( std::size_t k = (same_) ? 0 : cmdata::indexing::n_bins(cutoff_); k < cmdata::indexing::n_bins(cutoff_); k++ )
          {
            // sum+= dx_ * interm_same_maxcdf_mol_[i][ii][jj][k];
            // if (sum > 1.0) sum=1.0;
            // interm_same_maxcdf_mol_[i][ii][jj][k] = sum;
          }
          if (same_) 
          {
            for ( int jj = ii+1; jj < nresmol2_[i]; jj++ ) interm_same_mat_density_[i][jj][ii] = interm_same_mat_density_[i][ii][jj];
          }
          // if (same_) interm_same_maxcdf_mol_[i][jj][ii] = interm_same_maxcdf_mol_[i][ii][jj];
          // if (intra_) intram_mat_density_[i][jj][ii] = intram_mat_density_[i][ii][jj];
        // }
      }
      for (std::size_t j = i + 1; j < nresmol2_.size() && cross_; j++)
      {
        for (int ii = 0; ii < nresmol2_[i]; ii++)
        {
          for (int jj = 0; jj < nresmol2_[j]; jj++)
          {
            double inv_num_mol_cross = inv_num_mol_unique_[i];
            normalize_inter_cross(cross_index_[i][j], ii, norm, inv_num_mol_cross, interm_cross_mat_density_);
            // normalize_inter_cross(cross_index_[i][j], ii, jj, norm, inv_num_mol_cross, interm_cross_maxcdf_mol_);

            double sum = 0.0;
            for ( std::size_t k = (cross_) ? 0 : cmdata::indexing::n_bins(cutoff_); k < cmdata::indexing::n_bins(cutoff_); k++ )
            {
              // sum += dx_ * interm_cross_maxcdf_mol_[cross_index_[i][j]][ii][jj][k];
              // if (sum > 1.0) sum = 1.0;
              // interm_cross_maxcdf_mol_[cross_index_[i][j]][ii][jj][k] = sum;
            }
          }
        }
      }
    }
  }
  

static void f_write_inter_same(const std::string &output_prefix,
  std::size_t i, int ii, const std::vector<int> &natmol2,
  const std::vector<std::vector<std::vector<double>>> &interm_same_mat_density
)
{
  std::cout << "Writing data for atom " << ii << "..." << std::endl;
  std::filesystem::path ffh_inter = output_prefix + "inter_mol_" + std::to_string(i + 1) + "_" + std::to_string(i + 1) + "_aa_" + std::to_string(ii + 1) + ".dat";
  std::ofstream fp_inter(ffh_inter);
  for ( std::size_t k = 0; k < interm_same_mat_density[i][ii].size(); k++ )
  {
    fp_inter << " " << COUT_FLOAT_PREC6 << interm_same_mat_density[i][ii][k];
    fp_inter << "\n";
  }
  fp_inter.close();
}

static void f_write_inter_cross(const std::string &output_prefix,
  std::size_t i, std::size_t j, int ii, const std::vector<int> &natmol2,
  const std::vector<std::vector<int>> &cross_index,
  const cmdata_matrix &interm_cross_mat_density
)
{
  std::filesystem::path ffh = output_prefix + "rinter_mol_" + std::to_string(i + 1) + "_" + std::to_string(j + 1) + "_aa_" + std::to_string(ii + 1) + ".dat";
  std::ofstream fp(ffh);
  for ( std::size_t k = 0; k < interm_cross_mat_density[cross_index[i][j]][ii].size(); k++ )
  {
    for (int jj = 0; jj < natmol2[j]; jj++)
    {
      fp << " " << std::fixed << std::setprecision(7) << interm_cross_mat_density[cross_index[i][j]][ii][jj];
    }
    fp << "\n";
  }
  fp.close();
}

  void write_output( const std::string &output_prefix )
  {
    std::cout << "Writing data... " << std::endl;
    // using ftype_write_intra = cmdata::ftypes::function_traits<decltype(&cmdata::io::f_write_intra)>;
    using ftype_write_inter_same = cmdata::ftypes::function_traits<decltype(&f_write_inter_same)>;
    using ftype_write_inter_cross = cmdata::ftypes::function_traits<decltype(&f_write_inter_cross)>;
    std::function<ftype_write_inter_same::signature> write_inter_same = cmdata::ftypes::do_nothing<ftype_write_inter_same>();
    std::function<ftype_write_inter_cross::signature> write_inter_cross = cmdata::ftypes::do_nothing<ftype_write_inter_cross>();

    
    if (same_)
    {
      write_inter_same = f_write_inter_same;
      std::cout << "assos" << std::endl;
    }
    if (cross_) write_inter_cross = f_write_inter_cross;

    for (std::size_t i = 0; i < nresmol2_.size(); i++)
    {
      std::cout << "Writing data for molecule " << i << "..." << std::endl;
     
      std::filesystem::path ffh = output_prefix + "rinter_mol_" + std::to_string(i + 1) + "_" + std::to_string(i + 1) + ".dat";
      std::ofstream fp(ffh);

      cmdata::io::print_progress_bar(0.0);
      float progress = 0.0, new_progress = 0.0;    
      for (int ii = 0; ii < nresmol2_[i]; ii++)
      {
        new_progress = static_cast<float>(ii) / static_cast<float>(nresmol2_[i]);
        if (new_progress - progress > 0.01)
        {
          progress = new_progress;
          cmdata::io::print_progress_bar(progress);
        }
        // write_intra(output_prefix, i, ii, density_bins_, natmol2_, intram_mat_density_);
        // write_inter_same(output_prefix, i, ii, nresmol2_, interm_same_mat_density_);
        std::filesystem::path ffh = output_prefix + "rinter_mol_" + std::to_string(i + 1) + "_" + std::to_string(i + 1) + "_aa_" + std::to_string(ii + 1) + ".dat";
        
        for (int jj = 0; jj < nresmol2_[i]; jj++)
        {
          fp << ii << " " << jj << " " << std::fixed << std::setprecision(7) << interm_same_mat_density_[i][ii][jj] << "\n";
        }      
      }
      fp.close();
      for (std::size_t j = i + 1; j < nresmol2_.size(); j++)
      {
        for (int ii = 0; ii < nresmol2_[i]; ii++)
        {
          write_inter_cross(output_prefix, i, j, ii, nresmol2_, cross_index_, interm_cross_mat_density_);
        }
      }
    }
    cmdata::io::print_progress_bar(1.0);
    std::cout << "\nFinished!" << std::endl;
  }
};


} // namespace cmdata

#endif // _CM_DATA_HPP
