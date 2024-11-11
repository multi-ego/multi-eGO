#ifndef _RESDATA_RESDATA_HPP
#define _RESDATA_RESDATA_HPP

// gromacs includes
#include <gromacs/trajectoryanalysis/topologyinformation.h>
#include <gromacs/math/vec.h>
#include <gromacs/pbcutil/pbc.h>
#include <gromacs/fileio/tpxio.h>
#include <gromacs/fileio/confio.h>
#include <gromacs/topology/index.h>

// resdata includes
#include "io.hpp"
#include "indexing.hpp"
#include "parallel.hpp"
#include "density.hpp"
#include "block_avg.hpp"
//#include "mindist.hpp"
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
#include <chrono>


// xdrfile includes
#include <xdrfile.h>
#include <xdrfile_xtc.h>

namespace resdata
{

class RESData
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
  int nindex_pre_ndx;
  int nindex_=0;
  // std::vector<uint> selection_;
  gmx::RangePartitioning mols_pre_ndx;
  gmx::RangePartitioning mols_;
  std::vector<int> natmol2_;
  std::vector<int> mol_id_;
  std::vector<int> num_mol_unique_;
  bool use_index_ = false;
  // weights fields
  double weights_sum_;
  std::string weights_path_;
  std::vector<double> weights_;

  std::string traj_path_;
  // pbc fields
  bool no_pbc_;
  t_pbc *pbc_;
  PbcType pbcType_;

  // frame fields
  resdata::xtc::Frame *frame_;
  XDRFILE *trj_;

  //define the resmat matrix mol_i_type, mol_i, res_i, res_j for probability and distance
  using resmat_matrix = std::vector<std::vector<std::vector<double>>>;
  using resmat_matrix_temp = std::vector<std::vector<std::vector<std::vector<double>>>>;
  std::vector<int> residue_indeces_;
  resmat_matrix resmat_same_p_;
  resmat_matrix resmat_same_d_;
  resmat_matrix resmat_intra_p_;
  resmat_matrix resmat_intra_d_;
  resmat_matrix resmat_cross_p_;
  resmat_matrix resmat_cross_d_;
  std::vector<std::vector<int>> cross_index_;
  std::vector<int> num_res_per_molecule;

  resmat_matrix_temp appo_cross_d;
  resmat_matrix_temp appo_cross_p;
  resmat_matrix_temp appo_same_d;
  resmat_matrix_temp appo_same_p;  
  resmat_matrix_temp appo_intra_d;
  resmat_matrix_temp appo_intra_p;
  
  //BLOCK AVERAGE
  std::vector<int> block_sizes_;
  using block_resmat_matrix = std::vector<std::vector<std::vector<std::vector<double>>>>;
  block_resmat_matrix block_resmat_intra_p_;
  block_resmat_matrix block_resmat_intra_p2_;
  block_resmat_matrix block_resmat_intra_p_temp;
  block_resmat_matrix block_resmat_intra_std_;
  std::vector<std::vector<std::vector<std::vector<int>>>> block_counter_intra;
  block_resmat_matrix block_resmat_same_p_;
  block_resmat_matrix block_resmat_same_p2_;
  block_resmat_matrix block_resmat_same_p_temp;
  block_resmat_matrix block_resmat_same_std_;
  std::vector<std::vector<std::vector<std::vector<int>>>> block_counter_same;
  block_resmat_matrix block_resmat_cross_p_;
  block_resmat_matrix block_resmat_cross_p2_;
  block_resmat_matrix block_resmat_cross_p_temp;
  block_resmat_matrix block_resmat_cross_std_;
  std::vector<std::vector<std::vector<std::vector<int>>>> block_counter_cross;

  bool blk_avg_;

  t_topology top_ress_;
  rvec* x;
  int        natoms_res;
  int        isize_;
  int*       index_;

  // temporary containers for maxcdf operations
  std::vector<std::vector<std::mutex>> frame_same_mutex_;
  std::vector<std::vector<std::mutex>> frame_cross_mutex_;

  // parallelization fields
  int num_threads_;
  int num_mol_threads_;
  std::vector<std::thread> threads_;
  std::vector<std::thread> mol_threads_;
  resdata::parallel::Semaphore semaphore_;

  // mode selection, booleans and functions
  std::string mode_;
  bool intra_ = false, same_ = false, cross_ = false;

  static void molecule_routine(
    const int i, const int nindex_, t_pbc *pbc, rvec *x,
    const double cut_sig_2_, 
    const std::vector<int> &natmol2_, const std::vector<int> &num_mol_unique_, const std::vector<int> &mol_id_, 
    const std::vector<std::vector<int>> &cross_index_, 
    const double mcut2_, rvec *xcm_, const gmx::RangePartitioning &mols_, gmx_mtop_t *mtop_, 
    std::vector<std::vector<std::mutex>> &frame_cross_mutex_, 
    std::vector<std::vector<std::mutex>> &frame_same_mutex_, 
    resdata::parallel::Semaphore &semaphore_, double weight,const std::vector<int> &residue_indeces_, 
    resmat_matrix &resmat_same_d_, resmat_matrix &resmat_same_p_,
    resmat_matrix_temp &appo_same_d,    resmat_matrix_temp &appo_same_p,    
    resmat_matrix &resmat_intra_d_,resmat_matrix &resmat_intra_p_,
    resmat_matrix_temp &appo_intra_d,   resmat_matrix_temp &appo_intra_p,
    resmat_matrix &resmat_cross_d_,resmat_matrix &resmat_cross_p_,
    resmat_matrix_temp &appo_cross_d,   resmat_matrix_temp &appo_cross_p,
    const bool same, const bool intra, const bool cross,
    const int* index,
    const std::vector<int> &num_res_per_molecule
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
      if (tmp_i == num_mol_unique_.size()) 
      break;
    }
    if (mol_i == num_mol_unique_[mol_id_[i]]) mol_i = 0;
    
    int molb = 0;
    /* Loop over molecules  */
    for (int j = 0; j < nindex_; j++)
    {
      if((same || intra) && (!cross)){
        if(mol_id_[i]!=mol_id_[j]) continue;
      }
      if((!same && !intra) && cross){
        if(mol_id_[i]==mol_id_[j]){
          continue;
        }
      }
      if((!same && !intra && !cross)){
        continue;
      }
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
      std::size_t a_i = 0;
      // GMX_RELEASE_ASSERT(mols_.numBlocks() > 0, "Cannot access index[] from empty mols");

      /* cycle over the atoms of a molecule i */
      for (std::size_t iii = mols_.block(i).begin(); iii < mols_.block(i).end(); iii++)
      {
        int ii = index[static_cast<int>(iii)];
        std::size_t a_j = 0;
        mtopGetAtomAndResidueName(*mtop_, ii, &molb, &atomname, nullptr, nullptr, nullptr);    
        if (atomname[0] == 'H')
        {
          a_i++;
          continue;
        }
        /* cycle over the atoms of a molecule j */
        for (std::size_t jjj = mols_.block(j).begin(); jjj < mols_.block(j).end(); jjj++)
        {
          int jj = index[static_cast<int>(jjj)];
          mtopGetAtomAndResidueName(*mtop_, jj, &molb, &atomname, nullptr, nullptr, nullptr);
          if (atomname[0] == 'H')
          {
            a_j++;
            continue;
          }
          std::size_t delta  = iii-jjj-(i-j)*natmol2_[mol_id_[i]];
          rvec sym_dx;
          if (pbc != nullptr) pbc_dx(pbc, x[ii], x[jj], sym_dx);
          else rvec_sub(x[ii], x[jj], sym_dx);
          double dx2 = iprod(sym_dx, sym_dx);
          if (i==j) 
          {
            if(intra)
            {
              if (dx2 < cut_sig_2_)
              { // intra molecule species
                  std::size_t same_mutex_index = resdata::indexing::mutex_access(mol_id_[i], residue_indeces_[ii]-1, residue_indeces_[jj]-1, num_res_per_molecule);
                  std::unique_lock lock(frame_same_mutex_[mol_id_[i]][same_mutex_index]);
                  appo_intra_d[mol_id_[i]][mol_i][residue_indeces_[ii]-1][residue_indeces_[jj]-1] = std::min( appo_intra_d[mol_id_[i]][mol_i][residue_indeces_[ii]-1][residue_indeces_[jj]-1],std::sqrt(dx2));
                  appo_intra_p[mol_id_[i]][mol_i][residue_indeces_[ii]-1][residue_indeces_[jj]-1] = 1.;
                  lock.unlock();
              }
            }
          }
          else
          {
            if (mol_id_[i]==mol_id_[j])
            { // inter same molecule specie
              if(same)
              {
                if (dx2 < cut_sig_2_)
                {
                  std::size_t same_mutex_index = resdata::indexing::mutex_access(mol_id_[i], residue_indeces_[ii]-1, residue_indeces_[jj]-1, num_res_per_molecule);
                  std::unique_lock lock(frame_same_mutex_[mol_id_[i]][same_mutex_index]);
                  appo_same_d[mol_id_[i]][mol_i][residue_indeces_[ii]-1][residue_indeces_[jj]-1] = std::min( appo_same_d[mol_id_[i]][mol_i][residue_indeces_[ii]-1][residue_indeces_[jj]-1],std::sqrt(dx2));
                  appo_same_p[mol_id_[i]][mol_i][residue_indeces_[ii]-1][residue_indeces_[jj]-1] = 1.;
                  lock.unlock();
                }
                if (delta!=0.) {
                  // this is to account for inversion atom/molecule
                  int ii_delta =index[iii-delta];
                  int jj_delta =index[jjj+delta];

                  if (pbc != nullptr) pbc_dx(pbc, x[ii_delta], x[jj_delta], sym_dx);
                  else rvec_sub(x[ii_delta], x[jj_delta], sym_dx);
                  dx2 = iprod(sym_dx, sym_dx);
                  if (dx2 < cut_sig_2_)
                  {
                    std::size_t same_mutex_index = resdata::indexing::mutex_access(mol_id_[i],  residue_indeces_[ii]-1, residue_indeces_[jj]-1, num_res_per_molecule);
                    std::unique_lock lock(frame_same_mutex_[mol_id_[i]][same_mutex_index]);
                    appo_same_d[mol_id_[i]][mol_i][residue_indeces_[jj_delta]-1][residue_indeces_[ii_delta]-1] = std::min( appo_same_d[mol_id_[i]][mol_i][residue_indeces_[jj_delta]-1][residue_indeces_[ii_delta]-1],std::sqrt(dx2));
                    appo_same_p[mol_id_[i]][mol_i][residue_indeces_[jj_delta]-1][residue_indeces_[ii_delta]-1] = 1.;
                    lock.unlock();
                  }
                }
              }
            } 
            else
            { // inter cross molecule species
              if (dx2 < cut_sig_2_ && cross)
              {
                  std::size_t cross_mutex_index = resdata::indexing::mutex_access(mol_id_[j], residue_indeces_[ii]-1, residue_indeces_[jj]-1, num_res_per_molecule);
                  std::unique_lock lock(frame_cross_mutex_[cross_index_[mol_id_[i]][mol_id_[j]]][cross_mutex_index]);
                  appo_cross_d[cross_index_[mol_id_[i]][mol_id_[j]]][mol_i][residue_indeces_[ii]-1][residue_indeces_[jj]-1] = std::min( appo_cross_d[cross_index_[mol_id_[i]][mol_id_[j]]][mol_i][residue_indeces_[ii]-1][residue_indeces_[jj]-1],std::sqrt(dx2));
                  appo_cross_p[cross_index_[mol_id_[i]][mol_id_[j]]][mol_i][residue_indeces_[ii]-1][residue_indeces_[jj]-1] = 1.;
                  lock.unlock();
              }
            }
          }
          ++a_j;
        }
        ++a_i;
      }
      ++mol_j;
    }
    semaphore_.release();
  }

public:
  RESData(
    const std::string &top_path, const std::string &traj_path,
    double cutoff, double mol_cutoff, int nskip, int num_threads, int num_mol_threads,
    int dt, const std::string &mode, const std::string &weights_path, 
    bool no_pbc, float t_begin, float t_end, const std::string &index_path, bool blk_avg
  ) : cutoff_(cutoff), mol_cutoff_(mol_cutoff), nskip_(nskip), num_threads_(num_threads), num_mol_threads_(num_mol_threads),
      mode_(mode), weights_path_(weights_path),
      no_pbc_(no_pbc), dt_(dt), t_begin_(t_begin), t_end_(t_end), blk_avg_(blk_avg)
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

    traj_path_ = traj_path;
    int natom;
    long unsigned int nframe;
    int64_t *offsets;

    PbcType    pbcType;
    matrix            box = { { 0 } };
    //read_tps_conf(ftp2fn(efTPS, NFILE, fnm), &top_ress, &pbcType, &x, nullptr, box, FALSE);
    read_tps_conf(top_path.c_str(), &top_ress_, &pbcType, &x, nullptr, box, FALSE);

    frame_ = (resdata::xtc::Frame*)malloc(sizeof(resdata::xtc::Frame));
    std::cout << "Reading trajectory file " << traj_path << std::endl;
    read_xtc_header(traj_path.c_str(), &natom, &nframe, &offsets);
    *frame_ = resdata::xtc::Frame(natom);
    frame_->nframe = nframe;
    frame_->offsets = offsets;

    trj_ = xdrfile_open(traj_path.c_str(), "r");

    //define index
    if ( std::filesystem::exists(std::filesystem::path(index_path)) )
    {
      char*      grpname;
      use_index_ = true;
      get_index(&top_ress_.atoms, index_path.c_str(), 1, &isize_, &index_, &grpname);
    }
    else{
      index_ = new int[natom];
      for(int i=0; i<natom; i++){
        index_[i] = i;
        isize_ = natom;
      }
    }
    initAnalysis();
  }

  ~RESData()
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

    std::vector<int> natmol2_pre_ndx;
    std::vector<int> mol_id_pre_ndx;
    std::vector<int> num_mol_unique_pre_ndx;

    // get the number of atoms per molecule
    // equivalent to mols_ = gmx:gmx_mtop_molecules(*top.mtop());
    for (const gmx_molblock_t &molb : mtop_->molblock)
    {
      int natm_per_mol = mtop_->moltype[molb.type].atoms.nr;
      for (int i = 0; i < molb.nmol; i++) mols_pre_ndx.appendBlock(natm_per_mol);
      num_mol_unique_pre_ndx.push_back(molb.nmol);
    }

    // number of molecules
    nindex_pre_ndx = mols_pre_ndx.numBlocks();

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
    int n_mol_appo = 0, n_mol_counter=0,n_mol_true=0;
    int molb_index = 0;
    int natom_appo=0,natom_appo_2=0 ;

    //find the number of atom per molecule and the mol_id if there is no index
    for ( auto i : num_mol_unique_pre_ndx )
    {
      natmol2_pre_ndx.push_back(mols_pre_ndx.block(molb_index).end() - mols_pre_ndx.block(molb_index).begin());
      for ( int j = 0; j < i; j++ )
      {
        mol_id_pre_ndx.push_back(mol_id);
      }
      mol_id++;
      molb_index += i;
    }

    // find the true number of atoms and mol id per molecule if index is present
    printf("\nReading index file\n");
    mol_id = 0;
    molb_index = 0;

    for ( int i=0; i<nindex_pre_ndx ; i++ )
    {
      n_mol_appo = i-n_mol_counter;
      if(n_mol_appo==num_mol_unique_pre_ndx[mol_id] )
      {
        if(n_mol_true>0)num_mol_unique_.push_back(n_mol_true);

        n_mol_counter+=num_mol_unique_pre_ndx[mol_id] ;
        mol_id++;
        if(n_mol_true>0)
        {
          natmol2_.push_back(natom_appo_2);
        }
        n_mol_true = 0;
      }
      natom_appo=0;
      //check if all atoms in mols_.block(i) are in the index file
      for(int i_atom=mols_pre_ndx.block(i).begin(); i_atom<mols_pre_ndx.block(i).end(); i_atom++)
      {
        for(int k = 0; k<isize_;k++){
          if(index_[k]==i_atom)natom_appo++;
        }
      }
      if((natom_appo == 0))
      {
        printf("  X  mol: %i, mol id: %i, size pre ndx: %i, size after ndx: %i \n",i+1,mol_id, natmol2_pre_ndx[mol_id], natom_appo);
      }
      else if((natom_appo > 0) && (natom_appo != natmol2_pre_ndx[mol_id]))
      {
        printf("\nERROR: Number of atoms given by the index does not match the number of atoms of the corresponding molecule.\nIndex file should be used to remove hole molecules, otherwise do a post processing analysis\n\n");
        printf("Different number of atoms found in molecule %i (molecule type %i): number of atom selected in index: %i, number of atoms of molecule in topology: %i\n",i,mol_id, natom_appo, natmol2_pre_ndx[mol_id]);
        printf("Exit Code\n");
        exit(2);
      }
      else{
        printf("  V  mol: %i, mol id: %i, size pre ndx: %i, size after ndx: %i \n",i+1,mol_id, natmol2_pre_ndx[mol_id], natom_appo);
        n_mol_true +=1;
        nindex_++;
        mol_id_.push_back(mol_id);
        mols_.appendBlock(natom_appo);
        natom_appo_2 = natom_appo;
      }
      
      //molb_index += i;
    }
    if(n_mol_true>0){
      natmol2_.push_back(natom_appo_2);
      num_mol_unique_.push_back(n_mol_true);
    }

    // CREATE VECTOR of RESIDUE indeces to map atom to residue
    printf("\nAssigning residue index to atoms\n");
    DefineResidueIndexing(index_, mtop_, num_mol_unique_, natmol2_, num_mol_unique_pre_ndx, natmol2_pre_ndx, num_res_per_molecule, residue_indeces_);

    printf("\nNumber of different molecules %lu\n", natmol2_.size());
    bool check_same = false;
    for(std::size_t i=0; i<natmol2_.size();i++) {
      printf("    mol id: %lu, num mols: %u, size: %u,  nresidues:%u\n", i, num_mol_unique_[i], natmol2_[i],num_res_per_molecule[i]);
      if(num_mol_unique_[i]>1) check_same = true;
    }

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


    printf("\nActivating modes\n");
    if(!check_same && same_) same_ = false;
    if(cross_ && natmol2_.size()<2)
    {
      cross_ = false; 
      printf(":: deactivating cross mode (only 1 type of molecule found)\n");
    }

    if(nindex_>1 && (same_ || cross_))
    {
      printf("\n\n::::::::::::WARNING::::::::::::\nMore than 1 molcule found in the system.\nFix pbc before running resdata using pbc mol\n");
      printf(":::::::::::::::::::::::::::::::\n\n");
    }

    if (same_)
    {
      std::cout << ":: activating intermat same calculations" << std::endl;
      resmat_same_d_.resize(natmol2_.size());
      resmat_same_p_.resize(natmol2_.size());
      appo_same_d.resize(natmol2_.size());
      appo_same_p.resize(natmol2_.size());
      if(blk_avg_)
      {
      block_resmat_same_p_.resize(natmol2_.size() );
      block_resmat_same_p2_.resize(natmol2_.size());
      block_resmat_same_p_temp.resize(natmol2_.size() );
      block_resmat_same_std_.resize(natmol2_.size() );
      block_counter_same.resize(natmol2_.size());
      }
    }
    if (cross_)
    {
      std::cout << ":: activating intermat cross calculations" << std::endl;
      resmat_cross_d_.resize((natmol2_.size() * (natmol2_.size() - 1)) / 2);
      resmat_cross_p_.resize((natmol2_.size() * (natmol2_.size() - 1)) / 2);
      appo_cross_d.resize((natmol2_.size() * (natmol2_.size() - 1)) / 2);
      appo_cross_p.resize((natmol2_.size() * (natmol2_.size() - 1)) / 2);
      if(blk_avg_)
      {
      block_resmat_cross_p_.resize((natmol2_.size() * (natmol2_.size() - 1)) / 2);
      block_resmat_cross_p2_.resize((natmol2_.size() * (natmol2_.size() - 1)) / 2);
      block_resmat_cross_p_temp.resize((natmol2_.size() * (natmol2_.size() - 1)) / 2);
      block_resmat_cross_std_.resize((natmol2_.size() * (natmol2_.size() - 1)) / 2);
      block_counter_cross.resize((natmol2_.size() * (natmol2_.size() - 1)) / 2);
      }
    }
    if (intra_) 
    {
      std::cout << ":: activating intramat calculations" << std::endl;
      resmat_intra_d_.resize(natmol2_.size());
      resmat_intra_p_.resize(natmol2_.size());
      appo_intra_d.resize(natmol2_.size());
      appo_intra_p.resize(natmol2_.size());
      if(blk_avg_)
      {
      block_resmat_intra_p_.resize(natmol2_.size() );
      block_resmat_intra_p2_.resize(natmol2_.size());
      block_resmat_intra_p_temp.resize(natmol2_.size() );
      block_resmat_intra_std_.resize(natmol2_.size() );
      block_counter_intra.resize(natmol2_.size());
      }
    }

    if(blk_avg_)
    {
      int n_frames_blk = 0;
      int frnr_appo = 0;
      int natom_appo;
      long unsigned int nframe_appo;
      int64_t *offsets_appo;

      resdata::xtc::Frame *frame_appo;
      XDRFILE *trj_appo;

      frame_appo = (resdata::xtc::Frame*)malloc(sizeof(resdata::xtc::Frame));
      read_xtc_header(traj_path_.c_str(), &natom_appo, &nframe_appo, &offsets_appo);
      *frame_appo = resdata::xtc::Frame(natom_appo);
      frame_appo->nframe = nframe_appo;
      frame_appo->offsets = offsets_appo;
      trj_appo = xdrfile_open(traj_path_.c_str(), "r");
      printf("\nDefining Block sizes\n");

      while (frame_appo->read_next_frame(trj_appo, no_pbc_, pbcType_, pbc_) == exdrOK)
      {

        if ((frame_appo->time >= t_begin_ && (t_end_ < 0 || frame_appo->time <= t_end_ )) && // within time borders 
            ( dt_ == 0 || std::fmod(frame_appo->time, dt_) == 0) && (nskip_ == 0 || std::fmod(frnr_appo, nskip_) == 0)) // skip frames
        {
          n_frames_blk++;
        }
        frnr_appo++;
      }  
      
      printf("\nActivating Block average analysis with block sizes (in frames)\n");
      int i_blk = 4;
      while( i_blk<n_frames_blk/4 )
      {
        printf("%i  %i  ", i_blk, i_blk*3/2);
        block_sizes_.push_back(i_blk);
        block_sizes_.push_back(i_blk*3/2);
        i_blk *= 2;
      }
    printf("\n");
    }

    int cross_count = 0;
    if (cross_) cross_index_.resize(natmol2_.size(), std::vector<int>(natmol2_.size(), 0));
    for ( std::size_t i = 0; i < natmol2_.size(); i++ )
    {
      if (same_)
      {
        resmat_same_d_[i].resize(num_res_per_molecule[i], std::vector<double>(num_res_per_molecule[i], 0.));
        resmat_same_p_[i].resize(num_res_per_molecule[i], std::vector<double>(num_res_per_molecule[i], 0.));
        appo_same_d[i].resize(num_mol_unique_[i], std::vector<std::vector<double>>(num_res_per_molecule[i], std::vector<double>(num_res_per_molecule[i], 100.)));
        appo_same_p[i].resize(num_mol_unique_[i], std::vector<std::vector<double>>(num_res_per_molecule[i], std::vector<double>(num_res_per_molecule[i], 0.)));   
        if(blk_avg_)
        {
          block_resmat_same_p_[i].resize(block_sizes_.size(), std::vector<std::vector<double>>(num_res_per_molecule[i], std::vector<double>(num_res_per_molecule[i], 0.)));
          block_resmat_same_p2_[i].resize(block_sizes_.size(), std::vector<std::vector<double>>(num_res_per_molecule[i], std::vector<double>(num_res_per_molecule[i], 0.)));
          block_resmat_same_p_temp[i].resize(block_sizes_.size(), std::vector<std::vector<double>>(num_res_per_molecule[i], std::vector<double>(num_res_per_molecule[i], 0.)));
          block_resmat_same_std_[i].resize(block_sizes_.size(), std::vector<std::vector<double>>(num_res_per_molecule[i], std::vector<double>(num_res_per_molecule[i], 0.)));
          block_counter_same[i].resize(block_sizes_.size(), std::vector<std::vector<int>>(num_res_per_molecule[i], std::vector<int>(num_res_per_molecule[i], 0)));
        }   
      }
      if (intra_) 
      {
        resmat_intra_d_[i].resize(num_res_per_molecule[i], std::vector<double>(num_res_per_molecule[i], 0.));
        resmat_intra_p_[i].resize(num_res_per_molecule[i], std::vector<double>(num_res_per_molecule[i], 0.));
        appo_intra_d[i].resize(num_mol_unique_[i],std::vector<std::vector<double>>(num_res_per_molecule[i], std::vector<double>(num_res_per_molecule[i], 100.)));
        appo_intra_p[i].resize(num_mol_unique_[i], std::vector<std::vector<double>>(num_res_per_molecule[i], std::vector<double>(num_res_per_molecule[i], 0.)));
        if(blk_avg_)
        {
          block_resmat_intra_p_[i].resize(block_sizes_.size(), std::vector<std::vector<double>>(num_res_per_molecule[i], std::vector<double>(num_res_per_molecule[i], 0.)));
          block_resmat_intra_p2_[i].resize(block_sizes_.size(), std::vector<std::vector<double>>(num_res_per_molecule[i], std::vector<double>(num_res_per_molecule[i], 0.)));
          block_resmat_intra_p_temp[i].resize(block_sizes_.size(), std::vector<std::vector<double>>(num_res_per_molecule[i], std::vector<double>(num_res_per_molecule[i], 0.)));
          block_resmat_intra_std_[i].resize(block_sizes_.size(), std::vector<std::vector<double>>(num_res_per_molecule[i], std::vector<double>(num_res_per_molecule[i], 0.)));
          block_counter_intra[i].resize(block_sizes_.size(), std::vector<std::vector<int>>(num_res_per_molecule[i], std::vector<int>(num_res_per_molecule[i], 0)));
        }  
      }
      for ( std::size_t j = i + 1; j < natmol2_.size() && cross_; j++ )
      {
        resmat_cross_d_[cross_count].resize(num_res_per_molecule[i], std::vector<double>(num_res_per_molecule[j], 0.));
        resmat_cross_p_[cross_count].resize(num_res_per_molecule[i], std::vector<double>(num_res_per_molecule[j], 0.));
        appo_cross_d[cross_count].resize(num_mol_unique_[i],std::vector<std::vector<double>>(num_res_per_molecule[i], std::vector<double>(num_res_per_molecule[j], 0.)));
        appo_cross_p[cross_count].resize(num_mol_unique_[i], std::vector<std::vector<double>>(num_res_per_molecule[i], std::vector<double>(num_res_per_molecule[j], 0.)));
        
        if(blk_avg_)
        {
          block_resmat_cross_p_[cross_count].resize(block_sizes_.size(), std::vector<std::vector<double>>(num_res_per_molecule[i], std::vector<double>(num_res_per_molecule[j], 0.)));
          block_resmat_cross_p2_[cross_count].resize(block_sizes_.size(), std::vector<std::vector<double>>(num_res_per_molecule[i], std::vector<double>(num_res_per_molecule[j], 0.)));
          block_resmat_cross_p_temp[cross_count].resize(block_sizes_.size(), std::vector<std::vector<double>>(num_res_per_molecule[i], std::vector<double>(num_res_per_molecule[j], 0.)));
          block_resmat_cross_std_[cross_count].resize(block_sizes_.size(), std::vector<std::vector<double>>(num_res_per_molecule[i], std::vector<double>(num_res_per_molecule[j], 0.)));
          block_counter_cross[cross_count].resize(block_sizes_.size(), std::vector<std::vector<int>>(num_res_per_molecule[i], std::vector<int>(num_res_per_molecule[j], 0)));
        }
        cross_index_[i][j] = cross_count;
        cross_count++;
      }
    }

    mcut2_ = mol_cutoff_ * mol_cutoff_;
    cut_sig_2_ = (cutoff_) * (cutoff_ );
    xcm_ = (rvec*)malloc(nindex_ * sizeof(rvec));

    if (intra_ || same_) frame_same_mutex_.resize(natmol2_.size());
    if (cross_) frame_cross_mutex_.resize(cross_index_.size());
    std::size_t sum_cross_mol_sizes = 0;
    for ( std::size_t i = 0; i < natmol2_.size(); i++ )
    {
      if (intra_ || same_) frame_same_mutex_[i] = std::vector<std::mutex>(num_res_per_molecule[i] *  num_res_per_molecule[i]);
      for ( std::size_t j = i+1; j < natmol2_.size() && cross_; j++ )
      {
        frame_cross_mutex_[cross_index_[i][j]] = std::vector<std::mutex>(num_res_per_molecule[i] * num_res_per_molecule[j]);
      }
    }

    if (weights_path_ != "")
    {
      printf("Weights file provided. Reading weights from %s\n", weights_path_.c_str());
      weights_ = resdata::io::read_weights_file(weights_path_);
      printf("Found %li frame weights in file\n", weights_.size());
      double w_sum = std::accumulate(std::begin(weights_), std::end(weights_), 0.0, std::plus<>());
      printf("Sum of weights amounts to %lf\n", w_sum);
      weights_sum_ = 0.;
    }

    printf("Finished preprocessing. Starting frame-by-frame analysis.\n");
  }

  void run()
  {
    std::chrono::steady_clock::time_point begin;// = std::chrono::steady_clock::now();
    std::chrono::steady_clock::time_point end; //= std::chrono::steady_clock::now();
    std::cout << "Running frame-by-frame analysis" << std::endl;
    int frnr = 0;
    float progress = 0.0, new_progress = 0.0;
    resdata::io::print_progress_bar(progress);

    while (frame_->read_next_frame(trj_, no_pbc_, pbcType_, pbc_) == exdrOK)
    {
      new_progress = static_cast<float>(frnr) / static_cast<float>(frame_->nframe);
      if (new_progress - progress > 0.01)
      {
        progress = new_progress;
        resdata::io::print_progress_bar(progress);
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
        // #pragma omp parallel for num_threads(std::min(num_threads_, nindex_))
        for (int i = 0; i < nindex_; i++)
        {
          clear_rvec(xcm_[i]);
          double tm = 0.;
          for (int iii = mols_.block(i).begin(); iii < mols_.block(i).end(); iii++)
          {
            int ii = index_[iii];
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
        begin = std::chrono::steady_clock::now();
        for (int i = 0; i < nindex_; i++)
        {
          /* start molecule thread*/
          mol_threads_[i] = std::thread(
            molecule_routine, i, nindex_, pbc_, frame_->x, 
            cut_sig_2_, std::cref(natmol2_), std::cref(num_mol_unique_), std::cref(mol_id_), std::cref(cross_index_),
            mcut2_, xcm_, mols_, mtop_,
            std::ref(frame_cross_mutex_),std::ref(frame_same_mutex_),
            std::ref(semaphore_), weight, residue_indeces_, 
            std::ref(resmat_same_d_), std::ref(resmat_same_p_),std::ref(appo_same_d), std::ref(appo_same_p),
            std::ref(resmat_intra_d_), std::ref(resmat_intra_p_),std::ref(appo_intra_d), std::ref(appo_intra_p),
            std::ref(resmat_cross_d_), std::ref(resmat_cross_p_),std::ref(appo_cross_d), std::ref(appo_cross_p),
            same_, intra_, cross_, index_, std::cref(num_res_per_molecule)
          );
          /* end molecule thread */
        }
        for ( auto &thread : mol_threads_ ) thread.join();
        end = std::chrono::steady_clock::now();

        begin = std::chrono::steady_clock::now();
        for(int i = 0; i< natmol2_.size(); i++)
        { 
          if(intra_)
          {
            if(blk_avg_) resdata::blockavg::AccumulateIntraBlock(i,n_x_ ,std::cref(mol_id_), std::cref(natmol2_), std::ref(block_resmat_intra_p_), std::ref(block_resmat_intra_p2_), std::ref(block_resmat_intra_p_temp), std::ref(block_counter_intra),std::cref(block_sizes_), std::ref(appo_intra_p), std::cref(num_mol_unique_));
            resdata::density::AccumulateIntra(i, std::cref(mol_id_), std::ref(resmat_intra_d_), std::ref(resmat_intra_p_), std::ref(appo_intra_d), std::ref(appo_intra_p), std::cref(num_mol_unique_));
          }
          if(same_)
          {
            if(blk_avg_) resdata::blockavg::AccumulateInterSameBlock(i,n_x_ ,std::cref(mol_id_), std::cref(natmol2_), std::ref(block_resmat_same_p_), std::ref(block_resmat_same_p2_), std::ref(block_resmat_same_p_temp), std::ref(block_counter_same),std::cref(block_sizes_), std::ref(appo_same_p), std::cref(num_mol_unique_));
            resdata::density::AccumulateInterSame(i, std::cref(mol_id_), std::ref(resmat_same_d_), std::ref(resmat_same_p_), std::ref(appo_same_d), std::ref(appo_same_p), std::cref(num_mol_unique_));
          }
          if(cross_){
            if(blk_avg_) resdata::blockavg::AccumulateInterCrossBlock(i,n_x_ ,std::cref(mol_id_),std::cref(cross_index_), std::cref(natmol2_), std::ref(block_resmat_cross_p_), std::ref(block_resmat_cross_p2_), std::ref(block_resmat_cross_p_temp), std::ref(block_counter_cross),std::cref(block_sizes_), std::ref(appo_cross_p), std::cref(num_mol_unique_));
            resdata::density::AccumulateInterCross(i, std::cref(mol_id_),std::cref(cross_index_), std::cref(natmol2_), std::ref(resmat_cross_d_), std::ref(resmat_cross_p_), std::ref(appo_cross_d), std::ref(appo_cross_p),std::cref(num_mol_unique_));
          }
        }
        end = std::chrono::steady_clock::now();
        ++n_x_;
      }
      ++frnr;
    }

    resdata::io::print_progress_bar(1.0);
    std::cout<<"number of frames: "<<frnr<<std::endl;
    std::cout<<"frames analized: "<<n_x_<<std::endl;

    // NORMILIZE
    if(intra_)
    {
      resdata::density::NormilizeIntra(n_x_, std::ref(resmat_intra_d_) , std::ref(resmat_intra_p_), std::cref(num_mol_unique_));
      if(blk_avg_)resdata::blockavg::NormilizeIntraBlock(std::ref(block_resmat_intra_std_),std::ref(block_resmat_intra_p_),std::ref(block_resmat_intra_p2_), std::ref(block_counter_intra), std::cref(num_mol_unique_));
    }
    if(same_)
    {
      resdata::density::NormilizeInterSame(n_x_, std::ref(resmat_same_d_) , std::ref(resmat_same_p_), std::cref(num_mol_unique_));
      if(blk_avg_)resdata::blockavg::NormilizeInterSameBlock(std::ref(block_resmat_same_std_),std::ref(block_resmat_same_p_),std::ref(block_resmat_same_p2_), std::ref(block_counter_same), std::cref(num_mol_unique_));
    }
    if(cross_)
    {
      resdata::density::NormilizeInterCross(n_x_, std::ref(resmat_cross_d_) , std::ref(resmat_cross_p_));
      if(blk_avg_)resdata::blockavg::NormilizeInterCrossBlock(std::ref(block_resmat_cross_std_),std::ref(block_resmat_cross_p_),std::ref(block_resmat_cross_p2_), std::ref(block_counter_cross));
    }
  }

  void write_output( const std::string &output_prefix )
  {
    std::cout << "Writing data... " << std::endl;
    using ftype_write_intra_res = resdata::ftypes::function_traits<decltype(&resdata::io::f_write_intra_res)>;
    using ftype_write_inter_same_res = resdata::ftypes::function_traits<decltype(&resdata::io::f_write_inter_same_res)>;
    using ftype_write_inter_cross_res = resdata::ftypes::function_traits<decltype(&resdata::io::f_write_inter_cross_res)>;
    using ftype_write_inter_cross_res_block = resdata::ftypes::function_traits<decltype(&resdata::io::f_write_inter_cross_res_block)>;
    using ftype_write_inter_same_res_block = resdata::ftypes::function_traits<decltype(&resdata::io::f_write_inter_same_res_block)>;
    using ftype_write_intra_res_block = resdata::ftypes::function_traits<decltype(&resdata::io::f_write_intra_res_block)>;
    std::function<ftype_write_intra_res::signature> write_intra_res = resdata::ftypes::do_nothing<ftype_write_intra_res>();
    std::function<ftype_write_inter_same_res::signature> write_inter_same_res = resdata::ftypes::do_nothing<ftype_write_inter_same_res>();
    std::function<ftype_write_inter_cross_res::signature> write_inter_cross_res = resdata::ftypes::do_nothing<ftype_write_inter_cross_res>();
    std::function<ftype_write_inter_cross_res_block::signature> write_inter_cross_res_block = resdata::ftypes::do_nothing<ftype_write_inter_cross_res_block>();
    std::function<ftype_write_inter_same_res_block::signature> write_inter_same_res_block = resdata::ftypes::do_nothing<ftype_write_inter_same_res_block>();
    std::function<ftype_write_intra_res_block::signature> write_intra_res_block = resdata::ftypes::do_nothing<ftype_write_intra_res_block>();
    if (intra_) write_intra_res = resdata::io::f_write_intra_res;
    if (same_) write_inter_same_res = resdata::io::f_write_inter_same_res;
    if (cross_) write_inter_cross_res = resdata::io::f_write_inter_cross_res;
    if (cross_ && blk_avg_) write_inter_cross_res_block = resdata::io::f_write_inter_cross_res_block;
    if (same_ && blk_avg_) write_inter_same_res_block = resdata::io::f_write_inter_same_res_block;
    if (intra_ && blk_avg_) write_intra_res_block = resdata::io::f_write_intra_res_block;

    for (std::size_t i = 0; i < natmol2_.size(); i++)
    {
      std::cout << "Writing data for molecule " << i << "..." << std::endl;
      resdata::io::print_progress_bar(0.0);
      float progress = 0.0, new_progress = 0.0;    

      if(same_)
      {
        write_inter_same_res(output_prefix, i,  num_res_per_molecule,residue_indeces_, resmat_same_d_, resmat_same_p_);
        if(blk_avg_)write_inter_same_res_block(output_prefix, i,  num_res_per_molecule,residue_indeces_, resmat_same_p_, block_resmat_same_std_, block_sizes_);
      }
      if(intra_)
      {
        write_intra_res(output_prefix, i,  num_res_per_molecule, residue_indeces_, resmat_intra_d_, resmat_intra_p_);
        if(blk_avg_)write_intra_res_block(output_prefix, i,  num_res_per_molecule,residue_indeces_, resmat_intra_p_, block_resmat_intra_std_, block_sizes_);
      }
      for (std::size_t j = i + 1; j < natmol2_.size(); j++)
      {
         if(cross_) 
         {
          write_inter_cross_res(output_prefix, i, j,  num_res_per_molecule, cross_index_,residue_indeces_, resmat_cross_d_, resmat_cross_p_);
          if(blk_avg_)
          {
            printf("%i\n", block_sizes_.size());
            write_inter_cross_res_block(output_prefix, i, j,  num_res_per_molecule, cross_index_,residue_indeces_,resmat_cross_p_, block_resmat_cross_std_, block_sizes_);
          }
         }
      }
    }
    resdata::io::print_progress_bar(1.0);
    std::cout << "\nFinished!" << std::endl;
  }
};


} // namespace resdata

#endif // _RESDATA_HPP
