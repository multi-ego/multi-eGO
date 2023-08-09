/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2010- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * Implements gmx::analysismodules::CMData.
 *
 * \author multi-eGO development team
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "cmdata.h"

#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/trajectoryanalysis/topologyinformation.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/topology/mtop_lookup.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

#include <omp.h>
#include <thread>

#include <cmath>
#include <memory>
#include <string>
#include <algorithm>
#include <functional>
#include <numeric>
#include <sstream>
#include <fstream>

// #define timing
#ifdef timing
#include <chrono>
#define T_START(name) auto name##_start = std::chrono::high_resolution_clock::now()
static void _log_time(const int ms_time, const char* name ) { printf("timing :: %s :: %i us\n", name, ms_time);}
#define T_END(name) auto name##_end = std::chrono::high_resolution_clock::now(); auto name##_duration = std::chrono::duration_cast<std::chrono::microseconds>(name##_end - name##_start); _log_time(name##_duration.count(), #name)
#else
#define T_START(name) 0;
#define T_END(name) 0;
#endif // timing

// #define ACTIVATE_FILTER
#ifndef ACTIVATE_FILTER
#define FILTER(x) x
#else 
#define FILTER(x) 0
#endif


namespace gmx
{
 
namespace analysismodules
{
 
namespace
{

class CMData : public TrajectoryAnalysisModule
{
public:
  CMData();

  void initOptions(IOptionsContainer *options, TrajectoryAnalysisSettings *settings) override;
  void initAnalysis(const TrajectoryAnalysisSettings &settings, const TopologyInformation &top) override;

  void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc, TrajectoryAnalysisModuleData *pdata) override;

  void finishAnalysis(int nframes) override;
  void writeOutput() override;

private:
  Selection refsel_;
  bool histo_;
  double cutoff_;
  double mol_cutoff_;
  int n_x_;
  int nframe_;
  int max_mol_size_;
  int smax_mol_size_;
  rvec *xcm_ = nullptr;
  gmx_mtop_t *mtop_;
  std::vector<int> mol_id_;
  std::string outfile_inter_;
  std::string outfile_intra_;
  std::vector<double> inv_num_mol_;
  std::string sym_file_path_;
  std::vector<std::vector<std::vector<int>>> equivalence_list_;
  bool list_sym_;
  std::string list_sym_path_;

  std::vector<int> natmol2_;
  int nindex_;
  gmx::RangePartitioning mols_;
  std::vector<std::vector<int>> cross_index_;
  std::vector<t_atoms> molecules_;
  std::vector<double> density_bins_;
  int n_bins_;
  // std::vector<std::jthreads> threads_;
  int num_threads_;

  double mcut2_;
  double cut_sig_2_;
  std::vector<std::vector<std::vector<std::vector<double>>>> interm_same_mat_density_;
  std::vector<std::vector<std::vector<std::vector<double>>>> interm_cross_mat_density_;
  std::vector<std::vector<std::vector<std::vector<double>>>> intram_mat_density_;
  std::vector<std::vector<std::vector<std::vector<double>>>> interm_same_maxcdf_mol_;
  std::vector<std::vector<std::vector<std::vector<double>>>> interm_cross_maxcdf_mol_;

  // temporary containers for maxcdf operations
  std::vector<double> frame_same_mat_;
  std::vector<double> frame_cross_mat_;
};

CMData::CMData() : histo_(false),
                   outfile_intra_(""),
                   outfile_inter_(""),
                   sym_file_path_(""),
                   list_sym_(false) {}

void CMData::initOptions(IOptionsContainer *options, TrajectoryAnalysisSettings *settings)
{
  static const char *const desc[] = {
      "[THISMODULE] calculates the intra- and intermat properties for multi-eGO",
  };

  settings->setHelpText(desc);

  options->addOption(DoubleOption("cutoff")
                         .store(&cutoff_)
                         .defaultValue(0.75)
                         .description("Cutoff in which to consider contacts"));
  options->addOption(DoubleOption("mol_cutoff")
                         .store(&mol_cutoff_)
                         .defaultValue(6.0)
                         .description("Molecular cutoff in which to consider contacts intermolecularly"));
  options->addOption(FileNameOption("intra")
                         .store(&outfile_intra_)
                         .description("Output of the intra-molecular contacts"));
  options->addOption(FileNameOption("inter")
                         .store(&outfile_inter_)
                         .description("Output of the intra-molecular contacts"));
  options->addOption(BooleanOption("histo")
                         .store(&histo_)
                         .description("Set to true to output histograms"));
  options->addOption(SelectionOption("reference")
                         .store(&refsel_)
                         .required()
                         .dynamicMask()
                         .description("Groups to calculate distances to"));
  options->addOption(FileNameOption("sym")
                    .store(&sym_file_path_)
                    .description("Atoms symmetry file path"));
  options->addOption(BooleanOption("write_sym")
                         .store(&list_sym_)
                         .description("Write symmetry list"));

  // always require topology
  settings->setFlag(TrajectoryAnalysisSettings::efRequireTop);
}

static inline void read_symmetry_indices(
  const std::string &path,
  gmx_mtop_t &top,
  std::vector<std::vector<std::vector<int>>> &eq_list,
  const std::vector<int> &natmol2_,
  const std::vector<int> &start_index
  )
{

  eq_list.resize(natmol2_.size());
  for (std::size_t i = 0; i < natmol2_.size(); i++) {
    eq_list[i].resize(natmol2_[i]);
    for (int ii = 0; ii < natmol2_[i]; ii++)
    {
       eq_list[i][ii].push_back(ii);
    }
  }

  std::ifstream infile(path);
  if(path!=""&&!infile.good())
  {
    std::string errorMessage = "Cannot find the indicated symmetry file";
    GMX_THROW(InconsistentInputError(errorMessage.c_str()));
  }
  
  int molb = 0;
  std::string residue_entry, atom_entry_i, atom_entry_j;
  std::string line;
  // WARNING
  // this scales really bad... we should do the opposity and check for each atom pair in the same aminoacid if it has an equivalent atom
  while (std::getline(infile, line)) 
  {
    int atom1_index, atom2_index;
    std::istringstream iss(line);
    if (!(iss >> residue_entry >> atom_entry_i >> atom_entry_j)) // each necessary field is there
    {
      if (line=="") continue;
      printf("Skipping line\n%s\n due to syntax non-conformity\n", line.c_str());
      continue;
    }

    const char *atom_name_i, *atom_name_j, *residue_name_i, *residue_name_j;
    int a_i, a_j, resn_i, resn_j;
    for (std::size_t i = 0; i < natmol2_.size(); i++)
    {
      a_i = 0;
      for (int ii = start_index[i]; ii < start_index[i]+natmol2_[i]; ii++)
      {
        mtopGetAtomAndResidueName(top, ii, &molb, &atom_name_i, &resn_i, &residue_name_i, nullptr);
	      a_j = 0;
        for (int jj = start_index[i]; jj < start_index[i]+natmol2_[i]; jj++)
        {
          mtopGetAtomAndResidueName(top, jj, &molb, &atom_name_j, &resn_j, &residue_name_j, nullptr);
          if (((atom_name_i==atom_entry_i&&atom_name_j==atom_entry_j)||(atom_name_i==atom_entry_j&&atom_name_j==atom_entry_i))&&residue_entry==residue_name_i&&resn_i==resn_j)
          {
            bool insert = true;
            // check if element is already inserted
            for ( auto e : eq_list[i][a_i] )
            {
              if (e==a_j) insert = false;
            }
            // insert if not yet present
            if (insert) {
              eq_list[i][a_i].push_back(a_j);
            }
          }
	        a_j++;
        }
	      a_i++;
      }
    }
  }
}

// static inline void kernel_density_estimator(std::vector<double> &x, const std::vector<double> &bins, const double mu, const double norm)
// static inline void kernel_density_estimator(std::vector<double> &x, std::vector<double> &max_cdf, const std::vector<double> &bins, const double mu, const double norm)
static inline void kernel_density_estimator(std::vector<double>::iterator x, const std::vector<double> &bins, const double mu, const double norm)
{
  double h = 0.01;
  double from_x = std::max(mu - 2 * h, bins[0]);
  double to_x = std::min(mu + 2 * h, bins.back());
  auto is_geq_start = [&from_x](double i) { return i >= from_x; };
  auto is_geq_end = [&to_x](double i) { return i > to_x; };
  auto start = std::find_if(bins.begin(), bins.end(), is_geq_start);
  auto end = std::find_if(bins.begin(), bins.end(), is_geq_end);
  int from = std::distance(bins.begin(), start);
  int to = std::distance(bins.begin(), end);
  double scale = norm / (0.73853587 * h * std::sqrt(2. * M_PI));
  if (mu < h) scale *= 2.;
  double shift = std::exp(-2.);
  for (int i = from; i < to; i++)
  {
    double f = (mu - bins[i]) / h;
    double kernel = std::exp(-0.5 * f * f);
    x[i] += scale * (kernel - shift);
  }
}

static inline double calc_mean(const std::vector<double> &v, const double dx)
{
  double dm = 0.;
  double norm = 0.;
  for (auto it = v.begin(); it != v.end(); ++it)
  {
    unsigned i = std::distance(v.begin(), it);
    if (v[i] > 0.)
    {
      double d = (dx * static_cast<double>(i) + 0.5 * dx);
      dm += v[i] * d;
      norm += v[i];
    }
  }
  if (norm == 0.) norm = 1.;
  return dm / norm;
}

static inline double calc_prob(const std::vector<double> &v, const double dx)
{
  double prob = 0.;
  for (auto it = v.begin(); it != v.end(); ++it)
  {
    unsigned i = std::distance(v.begin(), it);
    if (v[i] > 0.) prob += v[i] * dx;
  }
  //if (prob > 1.) prob = 1.;
  return prob;
}

static inline int n_bins(const double cut, const double factor = 4.0)
{
  return cut / (0.01 / factor);
}

void CMData::initAnalysis(const TrajectoryAnalysisSettings &settings, const TopologyInformation &top)
{
  if (outfile_inter_ == "") outfile_inter_ = std::string("intermat.ndx");
  if (outfile_intra_ == "") outfile_intra_ = std::string("intramat.ndx");
  n_x_ = 0;
  nframe_ = 0;
  mtop_ = top.mtop();
  mols_ = gmx_mtop_molecules(*top.mtop());
  // number of molecules
  nindex_ = mols_.numBlocks();

  std::vector<int> num_mol;
  num_mol.push_back(1);
  int num_unique_molecules = 0;
  // number of atoms per molecule, assuming them identical when consecutive molecules have the same number of atoms
  natmol2_.push_back(mols_.block(0).end());
  for (int i = 1; i < nindex_; i++)
  {
    natmol2_.push_back(mols_.block(i).end() - mols_.block(i - 1).end());
    if (natmol2_[i] == natmol2_[i - 1]) num_mol[num_unique_molecules]++;
    else
    {
      num_mol.push_back(1);
      num_unique_molecules++;
    }
  }
  std::vector<int>::iterator it = std::unique(natmol2_.begin(), natmol2_.end());
  natmol2_.resize(std::distance(natmol2_.begin(), it));

  std::vector<int> start_index;
  mol_id_.push_back(0);
  start_index.push_back(0);
  num_unique_molecules = 0;
  inv_num_mol_.push_back(1. / (static_cast<double>(num_mol[num_unique_molecules])));

  for (int i = 1; i < nindex_; i++)
  {
    if (mols_.block(i).end() - mols_.block(i - 1).end() == natmol2_[num_unique_molecules])
    {
      start_index.push_back(start_index[i - 1]);
    }
    else
    {
      start_index.push_back(natmol2_[num_unique_molecules]);
      num_unique_molecules++;
    }
    mol_id_.push_back(num_unique_molecules);
    inv_num_mol_.push_back(1. / static_cast<double>(num_mol[num_unique_molecules]));
  }

  printf("number of different molecules %lu\n", natmol2_.size());
  for(std::size_t i=0; i<natmol2_.size();i++) printf("mol %lu num %u size %u\n", i, num_mol[i], natmol2_[i]);

  interm_same_mat_density_.resize(natmol2_.size());
  interm_cross_mat_density_.resize((natmol2_.size() * (natmol2_.size() - 1)) / 2);
  interm_same_maxcdf_mol_.resize(natmol2_.size());
  interm_cross_maxcdf_mol_.resize((natmol2_.size() * (natmol2_.size() - 1)) / 2);
  intram_mat_density_.resize(natmol2_.size());

  density_bins_.resize(n_bins(cutoff_));
  for (int i = 0; i < density_bins_.size(); i++)
    density_bins_[i] = cutoff_ / static_cast<double>(density_bins_.size()) * static_cast<double>(i) + cutoff_ / static_cast<double>(density_bins_.size() * 2);

  int cross_count = 0;
  cross_index_.resize(natmol2_.size(), std::vector<int>(natmol2_.size(), 0));
  for (std::size_t i = 0; i < natmol2_.size(); i++)
  {
    interm_same_mat_density_[i].resize(natmol2_[i], std::vector<std::vector<double>>(natmol2_[i], std::vector<double>(n_bins(cutoff_), 0)));
    interm_same_maxcdf_mol_[i].resize(natmol2_[i], std::vector<std::vector<double>>(natmol2_[i], std::vector<double>(n_bins(cutoff_), 0)));
    intram_mat_density_[i].resize(natmol2_[i], std::vector<std::vector<double>>(natmol2_[i], std::vector<double>(n_bins(cutoff_), 0)));
    for (std::size_t j = i + 1; j < natmol2_.size(); j++)
    {
      interm_cross_mat_density_[i].resize(natmol2_[i], std::vector<std::vector<double>>(natmol2_[j], std::vector<double>(n_bins(cutoff_), 0)));
      interm_cross_maxcdf_mol_[i].resize(natmol2_[i], std::vector<std::vector<double>>(natmol2_[j], std::vector<double>(n_bins(cutoff_), 0)));
      cross_index_[i][j] = cross_count;
      cross_count++;
    }
  }

  if (sym_file_path_=="") printf("No symmetry file provided. Running with standard settings.\n");
  else printf("Running with symmetry file %s\nReading file...\n", sym_file_path_.c_str());
  read_symmetry_indices(sym_file_path_, *mtop_, equivalence_list_, natmol2_, start_index);

  if (list_sym_)
  {
    printf("Writing out symmetry listing into %s\n", "sym_list.txt");
    std::fstream sym_list_file("sym_list.txt", std::fstream::out);
    for (int i = 0; i < equivalence_list_.size(); i++)
    {
      sym_list_file << "[ molecule_" << i << " ]\n";
      for (int j = 0; j < equivalence_list_[i].size(); j++)
      {
        sym_list_file << "atom " << j << ":";
        for (int k = 0; k < equivalence_list_[i][j].size(); k++)
        {
          sym_list_file << " " << equivalence_list_[i][j][k];
        }
        sym_list_file << "\n";
      }
      sym_list_file << "\n";
    }
    sym_list_file << "\n";
    sym_list_file.close();
  }

  n_bins_ = n_bins(cutoff_);

  mcut2_ = mol_cutoff_ * mol_cutoff_;
  cut_sig_2_ = (cutoff_ + 0.02) * (cutoff_ + 0.02);
  snew(xcm_, nindex_);

  int max_mol_size = *std::max_element(std::begin(natmol2_), std::end(natmol2_));
  int smax_mol_size = *std::max_element(std::begin(natmol2_), std::end(natmol2_),
                                        [max_mol_size](int a, int b) {
                                          if (a == max_mol_size) return true;
                                          if (b == max_mol_size) return false;                                         
                                          return a < b;
                                        });
  max_mol_size_ = max_mol_size;
  smax_mol_size_ = smax_mol_size;

  int sum_same_mol_sizes = 0; // std::accumulate(std::begin(natmol2_), std::end(natmol2_), 1., [](double acc, double v){ return acc += v * v;});
  int sum_cross_mol_sizes = 0;
  for (int i = 0; i < natmol2_.size(); i++)
  {
    sum_same_mol_sizes += natmol2_[i] * natmol2_[i];
    for (int j = i+1; j < natmol2_.size(); j++)
    {
      sum_cross_mol_sizes += natmol2_[i] * natmol2_[j];
    }
  }

  // frame_same_mat_.resize(natmol2_.size() * max_mol_size_ * max_mol_size_ * n_bins(cutoff_), 0);
  frame_same_mat_.resize(sum_same_mol_sizes * n_bins(cutoff_), 0);
  // frame_cross_mat_.resize(( ( natmol2_.size() * (natmol2_.size() - 1)) / 2 ) * max_mol_size_ * smax_mol_size_ * n_bins(cutoff_), 0);
  frame_cross_mat_.resize(sum_cross_mol_sizes * n_bins(cutoff_), 0);

  
  // num_threads_ = 12;
  // for (int i = 0; i < num_threads_; i++)
  // {
  //   threads_.emplace_back(

  //   );
  // }


  printf("Finished preprocessing. Starting frame-by-frame analysis.\n");
}

// #define access_same_(i, a_i, a_j) i * (natmol2_.size() * natmol2_[i] * natmol2_[i]) + a_i * (natmol2_[i] * natmol2_[i]) + a_j * natmol2_[i] 
#define access_same_(i, a_i, a_j) i * (natmol2_[i] * natmol2_[i] * n_bins_) + a_i * (natmol2_[i] * n_bins_) + a_j * n_bins_
#define access_cross_(i, j, a_i, a_j) cross_index_[i][j] * (natmol2_[i] * natmol2_[j] * n_bins_ ) + (natmol2_[j] * n_bins_) * a_i + n_bins_ * a_j
// #define access_cross_(i, j, a_i, a_j) cross_index_[i][j] * (( natmol2_.size() * (natmol2_.size() - 1)) / 2 ) * natmol2_[i] * natmol2_[j] + natmol2_[i] * natmol2_[j] * a_i + natmol2_[j] * a_j

static void accumulate_maxcdf_same(
  int start_im, const std::vector<int> start_i, const std::vector<int> start_j,
  int end_im, const std::vector<int> end_i, const std::vector<int> end_j,
  const int n_bins_, const std::vector<int> &natmol2_, const std::vector<double> &inv_num_mol,
  std::vector<double> &frame_same_mat,
  std::vector<std::vector<std::vector<std::vector<double>>>> &interm_same_mat_density,
  std::vector<std::vector<std::vector<std::vector<double>>>> &interm_same_maxcdf_mol
)
{
  bool first = true;
  int counter = 0;
  std::vector<bool> first_i(end_im - start_im, true);
  std::vector<bool> first_j(end_im - start_im, true);

  for (int im = start_im; im < end_im; im++)
  {
    FILTER(printf("OP ON im :: %i", im));
    int from_i = start_i[counter];
    int to_i = (im == end_im - 1) ? end_i[counter] : natmol2_[im];
    for (int i = from_i; i < to_i; i++)
    {
      int from_j = first_j[counter] ? start_j[counter] : i;
      int to_j = (i == end_i[counter]-1) ? end_j[counter] : natmol2_[im];
      for (int j = from_j; j < to_j; j++)
      {
        FILTER(printf("operating on (%i, %i, %i)\n", im, i, j));
        double sum = 0;
        int index = access_same_(im, i, j);
        for (int k = 0; k < n_bins_; ++k) 
        {
          sum+=frame_same_mat[index + k]*0.0025;
          if(sum>1.0) sum=1.0;
          interm_same_mat_density[im][i][j][k] += frame_same_mat[index + k]; 
          interm_same_maxcdf_mol[im][i][j][k] += sum*inv_num_mol[im]; 
        }
        interm_same_mat_density[im][j][i] = interm_same_mat_density[im][i][j];
        interm_same_maxcdf_mol[im][j][i] = interm_same_maxcdf_mol[im][i][j];
      }
      first_j[counter] = false;
    }
    first_i[counter] = false;
    counter++;
  }
  FILTER(printf("Thread finished\n"));
}

static void accumulate_maxcdf_cross(
  int start_im, const std::vector<int> start_jm, const std::vector<int> start_i, const std::vector<int> start_j,
  int end_im, const std::vector<int> end_jm, const std::vector<int> end_i, const std::vector<int> end_j,
  const int n_bins_, const std::vector<int> &natmol2_,
  const std::vector<std::vector<int>> &cross_index_,
  std::vector<double> &frame_cross_mat,
  std::vector<std::vector<std::vector<std::vector<double>>>> &interm_cross_mat_density,
  std::vector<std::vector<std::vector<std::vector<double>>>> &interm_cross_maxcdf_mol
)
{
  int jm_counter = 0;
  int counter = 0;
  // bool first_jm = true;
  // bool last_jm = false;
  // std::vector<bool> first_jm(start_jm.size(), true);
  // std::vector<bool> first_i(start_i.size(), true);
  // std::vector<bool> first_j(start_j.size(), true);

  for ( int im = start_im; im < end_im; im++ )
  {
    int jm_from = start_jm[jm_counter]; //(first_jm) ? start_jm[jm_counter] : (im + 1);
    int jm_to = end_jm[jm_counter]; // (last_jm) : end_jm[jm_counter] : end_jm[jm_counter];
    for ( int jm = jm_from; jm < natmol2_.size(); jm++ )
    {
      int from_i = start_i[counter];
      int to_i = end_i[counter];
      for (int i = from_i; i < to_i; i++)
      {
        int from_j = start_j[counter];
        int to_j = end_j[counter];
        for (int j = 0; j < natmol2_[j]; j++)
        {
          int index = access_cross_(im, jm, i, j);
          for (int kk = 0; kk < interm_cross_mat_density[cross_index_[im][jm]][i][0].size(); kk++) 
          {
            interm_cross_mat_density[cross_index_[im][jm]][i][j][kk] += frame_cross_mat[index + kk]; 
          }
        }
      }
      counter++;
    }
    // first_jm = false;
    jm_counter++;
  }
}

static void accumulate_max_cdf(
  const int n_bins_, const std::vector<int> &natmol2_, const std::vector<double> &inv_num_mol, // general parameters
  int start_im_same, const std::vector<int> start_i_same, const std::vector<int> start_j_same, // same parameters
  int end_im_same, const std::vector<int> end_i_same, const std::vector<int> end_j_same,
  std::vector<double> &frame_same_mat,
  std::vector<std::vector<std::vector<std::vector<double>>>> &interm_same_mat_density,
  std::vector<std::vector<std::vector<std::vector<double>>>> &interm_same_maxcdf_mol,
  int start_im_cross, const std::vector<int> start_jm_cross, const std::vector<int> start_i_cross, const std::vector<int> start_j_cross, // cross parameters
  int end_im_cross, const std::vector<int> end_jm_cross, const std::vector<int> end_i_cross, const std::vector<int> end_j_cross,
  const std::vector<std::vector<int>> &cross_index_,
  std::vector<double> &frame_cross_mat,
  std::vector<std::vector<std::vector<std::vector<double>>>> &interm_cross_mat_density,
  std::vector<std::vector<std::vector<std::vector<double>>>> &interm_cross_maxcdf_mol
)
{
  accumulate_maxcdf_same(
    start_im_same, start_i_same, start_j_same, end_im_same, end_i_same, end_j_same, n_bins_, natmol2_, inv_num_mol,
    frame_same_mat, interm_same_mat_density, interm_same_maxcdf_mol
  );
  // accumulate_maxcdf_cross(
  //   start_im_cross, start_jm_cross, start_i_cross, start_j_cross, end_im_cross, end_jm_cross, end_i_cross, end_j_cross,
  //   n_bins_, natmol2_, cross_index_, frame_cross_mat, interm_cross_mat_density, interm_cross_maxcdf_mol
  // );
}

void CMData::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc, TrajectoryAnalysisModuleData *pdata)
{
  // omp_set_nested(1);
  // WARNING IMPLEMENT
  int nskip = 0;
  // WARNING END

  // WARNING free memory again
  rvec *x = fr.x;

  if ((nskip == 0) || ((nskip > 0) && ((frnr % nskip) == 0)))
  {
    /* calculate the center of each molecule */
    for (int i = 0; (i < nindex_); i++)
    {
      clear_rvec(xcm_[i]);
      double tm = 0.;
      for (int ii = mols_.block(i).begin(); ii < mols_.block(i).end(); ii++)
      {
        for (int m = 0; (m < DIM); m++)
        {
          xcm_[i][m] += x[ii][m];
        }
        tm += 1.0;
      }
      for (int m = 0; (m < DIM); m++)
      {
        xcm_[i][m] /= tm;
      }
    }

    const char * atomname;

    /* Loop over molecules */
    for (int i = 0; i < nindex_; i++)
    {
      FILTER(printf("RUNNING %i\n", i));
      #ifdef timing
      auto start_zeroing = std::chrono::high_resolution_clock::now();
      #endif
      // works at time of 3500 ms for ttr
      // std::fill(std::begin(frame_same_mat_), std::end(frame_same_mat_), 0.);
      // std::fill(std::begin(frame_cross_mat_), std::end(frame_cross_mat_), 0.);

      #pragma omp parallel for num_threads(4)
      for ( int n = 0; n < frame_same_mat_.size(); n++ ) frame_same_mat_[n] = 0.;
      #pragma omp parallel for num_threads(4)
      for ( int n = 0; n < frame_cross_mat_.size(); n++ ) frame_cross_mat_[n] = 0.;

      #ifdef timing
      auto end_zeroing = std::chrono::high_resolution_clock::now();
      auto duration_zeroing = std::chrono::duration_cast<std::chrono::microseconds>(end_zeroing - start_zeroing);
      printf("frame :: %i, allocation took %li ms\n", frnr, duration_zeroing.count());
      #endif // timing
      
      // printf("acc %f\n", std::accumulate(std::begin(frame_same_mat_), std::end(frame_same_mat_), 0.));

      int molb = 0;
      /* Loop over molecules  */
      for (int j = 0; j < nindex_; j++)
      {
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

        int a_i = 0;
        GMX_RELEASE_ASSERT(mols_.numBlocks() > 0, "Cannot access index[] from empty mols");
        /* cycle over the atoms of a molecule i */
        for (int ii = mols_.block(i).begin(); ii < mols_.block(i).end(); ii++)
        {
          int a_j = 0;
          mtopGetAtomAndResidueName(*mtop_, ii, &molb, &atomname, nullptr, nullptr, nullptr);
          if (atomname[0] == 'H')
          {
            a_i++;
            continue;
          }
          /* cycle over the atoms of a molecule j */
          for (int jj = mols_.block(j).begin(); jj < mols_.block(j).end(); jj++)
          {
            mtopGetAtomAndResidueName(*mtop_, jj, &molb, &atomname, nullptr, nullptr, nullptr);
            if (atomname[0] == 'H')
            {
              a_j++;
              continue;
            }
            // check for chemical equivalence
            for (int eq_i = 0; eq_i < equivalence_list_[mol_id_[i]][a_i].size(); eq_i++)
            {
              for (int eq_j = 0; eq_j < equivalence_list_[mol_id_[j]][a_j].size(); eq_j++)
              {
                // get molecule-wise atom index considering equivalence
                int eqa_i  = equivalence_list_[mol_id_[i]][a_i][eq_i];             // molecule-wise equivalence index i
                int geqa_i = ii + (eqa_i - equivalence_list_[mol_id_[i]][a_i][0]); // global equivalence index i
                int eqa_j  = equivalence_list_[mol_id_[j]][a_j][eq_j];             // molecule-wise equivalence index j
                int geqa_j = jj + (eqa_j - equivalence_list_[mol_id_[j]][a_j][0]); // global equivalence index j
                int delta  = eqa_i - eqa_j;
                double nsym = static_cast<double>(equivalence_list_[mol_id_[i]][a_i].size()*equivalence_list_[mol_id_[i]][a_j].size());
                rvec sym_dx;
                if (pbc != nullptr) pbc_dx(pbc, x[geqa_i], x[geqa_j], sym_dx);
                else rvec_sub(x[geqa_i], x[geqa_j], sym_dx);
                double dx2 = iprod(sym_dx, sym_dx);
                if(i==j) 
                {
                  if (dx2 < cut_sig_2_)
                  {
                    #ifdef timing
                    auto start_intra_kde = std::chrono::high_resolution_clock::now();
                    #endif
                    kernel_density_estimator(std::begin(intram_mat_density_[mol_id_[i]][a_i][a_j]), density_bins_, std::sqrt(dx2), inv_num_mol_[i]/nsym);
                    #ifdef timing
                    auto end_intra_kde = std::chrono::high_resolution_clock::now();
                    auto duration_intra_kde = std::chrono::duration_cast<std::chrono::microseconds>(end_intra_kde - start_intra_kde);
                    printf("frame :: %i, intra took %li ms\n", frnr, duration_intra_kde.count());
                    #endif
                  }
                }
                else
                {
                  if(mol_id_[i]==mol_id_[j])
                  { // inter same molecule specie
                    if (dx2 < cut_sig_2_)
                    {
                      #ifdef timing
                      auto start_inter_kde = std::chrono::high_resolution_clock::now();
                      #endif
                      // kernel_density_estimator(frame_same_mat_[mol_id_[i]][a_i][a_j], density_bins_, std::sqrt(dx2), 1.0);
                      std::vector<double>::iterator starting_point = std::begin(frame_same_mat_) + access_same_(mol_id_[i], a_i, a_j);
                      kernel_density_estimator(starting_point, density_bins_, std::sqrt(dx2), 1.0);
                      #ifdef timing
                      auto end_inter_kde = std::chrono::high_resolution_clock::now();
                      auto duration_inter_kde = std::chrono::duration_cast<std::chrono::microseconds>(end_inter_kde - start_inter_kde);
                      printf("frame :: %i, inter took %li ms\n", frnr, duration_inter_kde.count());
                      #endif
                    }
                    if(delta!=0.) {
                      // this is to account for inversion atom/molecule
                      if (pbc != nullptr) pbc_dx(pbc, x[geqa_i-delta], x[geqa_j+delta], sym_dx);
                      else rvec_sub(x[geqa_i-delta], x[geqa_j+delta], sym_dx);
                      dx2 = iprod(sym_dx, sym_dx);
                      if (dx2 < cut_sig_2_)
                      {
                        #ifdef timing
                        auto start_interinv_kde = std::chrono::high_resolution_clock::now();
                        #endif
                        std::vector<double>::iterator starting_point = std::begin(frame_same_mat_) + access_same_(mol_id_[i], a_i, a_j);
                        kernel_density_estimator(starting_point, density_bins_, std::sqrt(dx2), 1.0);
                        #ifdef timing
                        auto end_interinv_kde = std::chrono::high_resolution_clock::now();
                        auto duration_interinv_kde = std::chrono::duration_cast<std::chrono::microseconds>(end_interinv_kde - start_interinv_kde);
                        printf("frame :: %i, interinv took %li ms\n", frnr, duration_interinv_kde.count());
                        #endif
                      }
                    }
                  } 
                  else
                  { // inter cross molecule specie
                    if (dx2 < cut_sig_2_)
                    {
                      #ifdef timing
                      auto start_cross_kde = std::chrono::high_resolution_clock::now();
                      #endif
                      std::vector<double>::iterator starting_point = std::begin(frame_cross_mat_) + access_cross_(mol_id_[i], mol_id_[j], a_i, a_j);
                      kernel_density_estimator(starting_point, density_bins_, std::sqrt(dx2), std::max(inv_num_mol_[i],inv_num_mol_[j])/nsym);
                      // frame_cross_count_[access_cross_(i, j, a_i, a_j)]++;
                      #ifdef timing
                      auto end_cross_kde = std::chrono::high_resolution_clock::now();
                      auto duration_cross_kde = std::chrono::duration_cast<std::chrono::microseconds>(end_cross_kde - start_cross_kde);
                      printf("frame :: %i, cross took %li ms\n", frnr, duration_cross_kde.count());
                      #endif
                    }
                  }
                }
              }
            }
            a_j++;
          }
          a_i++;
        }
      }

      /* accumulate the mean saturated cdf per molecule */
      ////////////////////////////////////////////////////
      //      TODO make it one continuous thing         //
      ////////////////////////////////////////////////////
      FILTER(printf("Accumulating run for molecule %i\n", i));
      int num_threads = 2;
      std::vector<std::thread> threads(num_threads);

      int num_ops_same = 0;
      for (int im = 0; im < natmol2_.size(); im++ ) num_ops_same += ( natmol2_[im] * ( natmol2_[im] + 1 ) ) / 2;
      int n_per_thread_same = num_ops_same / num_threads;
      int n_threads_same_uneven = num_ops_same % num_threads;
      int start_im_same = 0, end_im_same = 1; 
      std::vector<int> start_i_same({0}), start_j_same({0}), end_i_same({0}), end_j_same({0});
      int num_ops_cross = 0;
      for ( int im = 0; im < natmol2_.size(); im++ )
      {
        for ( int jm = im + 1; jm < natmol2_.size(); jm++ )
        {
          num_ops_cross += natmol2_[im] * natmol2_[jm];
        }
      }
      int n_per_thread_cross = num_ops_cross / num_threads;
      int n_threads_cross_uneven = num_ops_cross % num_threads;

      // std::vector<std::thread> threads_cross;
      int start_im_cross = 0, end_im_cross = 1;
      std::vector<int> start_jm_cross({start_im_cross + 1}), end_jm_cross({start_im_cross + 2}), start_i_cross({0}), end_i_cross({0}), start_j_cross({0}), end_j_cross({0});
      for ( int tid = 0; tid < num_threads; tid++ )
      {
        FILTER(printf("starting calculation for thread id :: %i of %i ", tid, num_threads));
        /* calculate same indices */
        int n_loop_operations_same = n_per_thread_same + (tid < n_threads_same_uneven ? 1 : 0);
        while (n_loop_operations_same - natmol2_[end_im_same - 1] + end_j_same.back() >= 0)
        {
          FILTER(printf("Entering while\n"));
          int sub_same = natmol2_[end_im_same - 1] - end_j_same.back();
          // FILTER(printf("sub :: %i\n", sub_same));
          n_loop_operations_same -= sub_same;
          end_i_same.back()++;
          end_j_same.back() = end_i_same.back();
          // FILTER(printf("im = (%i, %i), i = (%i, %i), j = (%i, %i), sub = %i, n_r = %i\n", 
          // start_im_same, end_im_same, start_i_same.back(), end_i_same.back(), start_j_same.back(), end_j_same.back(), sub, n_loop_operations));
          if ( end_i_same.back() == natmol2_[end_im_same - 1] )
          {
            FILTER(printf("increment i to %i and incresing vector sizes\n", start_im_same+1));
            end_im_same++;
            start_i_same.push_back(0);
            start_j_same.push_back(0);
            end_i_same.push_back(0);
            end_j_same.push_back(0);
          }
          if (n_loop_operations_same == 0) break;
        }
        end_j_same.back() += n_loop_operations_same;  
        FILTER(printf("Final coordinates im = (%i, %i), i = (%i, %i) & j = (%i, %i)\n", start_im_same, end_im_same, start_i_same.front(), end_i_same.back(), start_j_same.front(), end_j_same.back()));
        /* calculate cross indices */
        int n_loop_operations_cross = n_per_thread_cross + (tid < n_threads_cross_uneven ? 1 : 0);
        while ( n_loop_operations_cross - natmol2_[end_jm_cross.back() - 1] + end_j_cross.back() >= 0 )
        {
          FILTER(printf("entering cross while\n"));
          int sub_cross = natmol2_[end_jm_cross.back() - 1] - end_j_cross.back();
          n_loop_operations_cross -= sub_cross;
          end_i_cross.back()++;
          end_j_cross.back() = end_i_cross.back();
          if ( end_i_cross.back() == natmol2_[end_im_cross - 1] )
          {
            FILTER(printf("increasing cross vector sizes\n"));
            end_jm_cross.back()++;
            start_i_cross.push_back(0);
            start_j_cross.push_back(0);
            end_i_cross.push_back(0);
            end_j_cross.push_back(0);
          }
          if ( end_jm_cross.back() == natmol2_.size() )
          {
            FILTER(printf("increasing cross vector sizes 2\n"));
            end_im_cross++;
            start_jm_cross.push_back(end_im_cross);
            end_jm_cross.push_back(end_im_cross+1);
          }
          if (n_loop_operations_cross == 0) break;
        }
        end_j_cross.back() += n_loop_operations_cross;
        /* start thread same */
        // threads_same[tid] = std::thread(
        //   accumulate_maxcdf_same, start_im_same, start_i, start_j,
        //   end_im_same, end_i_same, end_j_same, n_bins_,
        //   std::cref(natmol2_), std::cref(inv_num_mol_),
        //   std::ref(frame_same_mat_), std::ref(interm_same_mat_density_),
        //   std::ref(interm_same_maxcdf_mol_)
        // );
        /* start thread */
        threads[tid] = std::thread(
          accumulate_max_cdf, n_bins_,  std::cref(natmol2_), std::cref(inv_num_mol_), // general parameters
          start_im_same, start_i_same, start_j_same, end_im_same, end_i_same, end_j_same, // same parameters
          std::ref(frame_same_mat_), std::ref(interm_same_mat_density_), std::ref(interm_same_maxcdf_mol_),
          start_im_cross, start_jm_cross, start_i_cross, start_j_cross, end_im_cross, end_jm_cross, // cross parameters 
          end_i_cross, end_j_cross, std::cref(cross_index_),
          std::ref(frame_cross_mat_),std::ref(interm_cross_mat_density_), std::ref(interm_cross_maxcdf_mol_)
        );
        /* end thread */
        // threads[tid].join();
        FILTER(printf("resetting the start and end values\n"));
        /* set new starts */
        start_im_same = end_im_same-1;
        start_i_same = std::vector<int>({end_i_same.back() - 1});
        start_j_same = std::vector<int>({end_j_same.back()});
        end_i_same = std::vector<int>({end_i_same.back()});
        end_j_same = std::vector<int>({end_j_same.back()});

        start_im_cross = end_im_cross - 1;
        start_jm_cross = std::vector<int>({end_jm_cross.back() - 1});
        start_i_cross = std::vector<int>({end_i_cross.back() - 1});
        start_j_cross = std::vector<int>({end_j_cross.back()});
        end_jm_cross = std::vector<int>({end_jm_cross.back()});
        end_i_cross = std::vector<int>({end_i_cross.back()});
        end_j_cross = std::vector<int>({end_j_cross.back()});
      }
      FILTER(printf("joining same threads\n"));
      // for ( auto &thread : threads_same ) thread.join();
      for ( auto &thread : threads ) thread.join();
      /* cross inter */
      // FILTER(printf("cross operation\n"));
      // int num_ops_cross = 0;
      // // int num_threads = 1;
      // for ( int im = 0; im < natmol2_.size(); im++ )
      // {
      //   for ( int jm = im + 1; jm < natmol2_.size(); jm++ )
      //   {
      //     num_ops_cross += natmol2_[im] * natmol2_[jm];
      //   }
      // }
      // int n_per_thread_cross = num_ops_cross / num_threads;
      // int n_threads_cross_uneven = num_ops_cross % num_threads;

      // std::vector<std::thread> threads_cross;
      // int start_im_cross = 0, end_im_cross = 1;
      // std::vector<int> start_jm_cross({start_im_cross + 1}), end_jm_cross({start_im_cross + 2}), start_i_cross({0}), end_i_cross({0}), start_j_cross({0}), end_j_cross({0});
      // for ( int tid = 0; tid < num_threads; tid = 0 )
      // {
      //   int n_loop_operations_cross = n_per_thread_cross + (tid < n_threads_cross_uneven ? 1 : 0);
      //   while ( n_loop_operations_cross - natmol2_[end_jm_cross.back() - 1] + end_j_cross.back() >= 0 )
      //   {
      //     int sub_cross = natmol2_[end_jm_cross.back() - 1] - end_j_cross.back();
      //     n_loop_operations_cross -= sub_cross;
      //     end_i_cross.back()++;
      //     end_j_cross.back() = end_i_cross.back();
      //     if ( end_i_cross.back() == natmol2_[end_im_cross - 1] )
      //     {
      //       end_jm_cross.back()++;
      //       start_i_cross.push_back(0);
      //       start_j_cross.push_back(0);
      //       end_i_cross.push_back(0);
      //       end_j_cross.push_back(0);
      //     }
      //     if ( end_jm_cross.back() == natmol2_.size() )
      //     {
      //       end_im_cross++;
      //       start_jm_cross.push_back(end_im_cross);
      //       end_jm_cross.push_back(end_im_cross+1);
      //     }
      //     if (n_loop_operations_cross == 0) break;
      //   }
      //   end_j_cross.back() += n_loop_operations_cross;

      //   /* start thread cross */
      //   threads_cross.emplace_back(
      //     std::thread(
      //       accumulate_maxcdf_cross, start_im_cross, start_jm_cross, start_i_cross, start_j_cross,
      //       end_im_cross, end_jm_cross, end_i_cross, end_j_cross, n_bins_,
      //       std::cref(natmol2_),//std::cref(inv_num_mol_),
      //       std::cref(cross_index_),
      //       std::ref(frame_cross_mat_), std::ref(interm_cross_mat_density_),
      //       std::ref(interm_cross_maxcdf_mol_)
      //     )
      //   );
      //   /* end thread cross */
      //   start_im_cross = end_im_cross - 1;
      //   start_jm_cross = std::vector<int>(end_jm_cross.back() - 1);
      //   start_i_cross = std::vector<int>(end_i_cross.back() - 1);
      //   start_j_cross = std::vector<int>(end_j_cross.back());
      //   end_jm_cross = std::vector<int>(end_jm_cross.back());
      //   end_i_cross = std::vector<int>(end_i_cross.back());
      //   end_j_cross = std::vector<int>(end_j_cross.back());
      // }
      // for ( auto &thread : threads_cross ) thread.join();
    }
    n_x_++;
  }
  nframe_++;
}

void CMData::finishAnalysis(int /*nframes*/)
{
  // normalisations
  double norm = 1. / n_x_;

  for (std::size_t i = 0; i < natmol2_.size(); i++)
  {
    for (int ii = 0; ii < natmol2_[i]; ii++)
    {
      for (int jj = ii; jj < natmol2_[i]; jj++)
      {
        std::transform(interm_same_maxcdf_mol_[i][ii][jj].begin(), 
                       interm_same_maxcdf_mol_[i][ii][jj].end(), 
                       interm_same_maxcdf_mol_[i][ii][jj].begin(), 
                       [&norm](auto &c) { return c * norm; });
        std::transform(interm_same_mat_density_[i][ii][jj].begin(), 
                       interm_same_mat_density_[i][ii][jj].end(), 
                       interm_same_mat_density_[i][ii][jj].begin(), 
                       [&norm](auto &c) { return c * norm; });
        std::transform(intram_mat_density_[i][ii][jj].begin(), 
                       intram_mat_density_[i][ii][jj].end(),
                       intram_mat_density_[i][ii][jj].begin(),
                       [&norm](auto &c) { return c * norm; });
        interm_same_mat_density_[i][jj][ii] = interm_same_mat_density_[i][ii][jj];
        interm_same_maxcdf_mol_[i][jj][ii] = interm_same_maxcdf_mol_[i][ii][jj];
        intram_mat_density_[i][jj][ii] = intram_mat_density_[i][ii][jj];
      }
    }
    for (std::size_t j = i + 1; j < natmol2_.size(); j++)
    {
      for (int ii = 0; ii < natmol2_[i]; ii++)
      {
        for (int jj = 0; jj < natmol2_[j]; jj++)
        {
          std::transform(interm_cross_mat_density_[cross_index_[i][j]][ii][jj].begin(),
                         interm_cross_mat_density_[cross_index_[i][j]][ii][jj].end(),
                         interm_cross_mat_density_[cross_index_[i][j]][ii][jj].begin(),
                         [&norm](auto &c) { return c * norm; });
        }
      }
    }
  }
}

void CMData::writeOutput()
{
  if (histo_)
  {
    for (std::size_t i = 0; i < natmol2_.size(); i++)
    {
      for (int ii = 0; ii < natmol2_[i]; ii++)
      {
        FILE *fp_inter = nullptr;
        FILE *fp_inter_cum = nullptr;
        FILE *fp_intra = nullptr;
        std::string ffh_inter = "inter_mol_" + std::to_string(i + 1) + "_" + std::to_string(i + 1) + "_aa_" + std::to_string(ii + 1) + ".dat";
        fp_inter = gmx_ffopen(ffh_inter, "w");
        std::string ffh_inter_cum = "inter_mol_c_" + std::to_string(i + 1) + "_" + std::to_string(i + 1) + "_aa_" + std::to_string(ii + 1) + ".dat";
        fp_inter_cum = gmx_ffopen(ffh_inter_cum, "w");
        std::string ffh_intra = "intra_mol_" + std::to_string(i + 1) + "_" + std::to_string(i + 1) + "_aa_" + std::to_string(ii + 1) + ".dat";
        fp_intra = gmx_ffopen(ffh_intra, "w");
        for (int k = 0; k < interm_same_mat_density_[i][ii][0].size(); k++)
        {
          fprintf(fp_inter, "%lf", density_bins_[k]);
          fprintf(fp_inter_cum, "%lf", density_bins_[k]);
          fprintf(fp_intra, "%lf", density_bins_[k]);
          for (int jj = 0; jj < natmol2_[i]; jj++)
          {
            fprintf(fp_inter, " %lf", interm_same_mat_density_[i][ii][jj][k]);
            fprintf(fp_inter_cum, " %lf", interm_same_maxcdf_mol_[i][ii][jj][k]);
            fprintf(fp_intra, " %lf", intram_mat_density_[i][ii][jj][k]);
          }
          fprintf(fp_inter, "\n");
          fprintf(fp_inter_cum, "\n");
          fprintf(fp_intra, "\n");
        }
        gmx_ffclose(fp_inter);
        gmx_ffclose(fp_inter_cum);
        gmx_ffclose(fp_intra);
      }
      for (std::size_t j = i + 1; j < natmol2_.size(); j++)
      {
        for (int ii = 0; ii < natmol2_[i]; ii++)
        {
          FILE *fp = nullptr;
          std::string ffh = "inter_mol_" + std::to_string(i + 1) + "_" + std::to_string(j + 1) + "_aa_" + std::to_string(ii + 1) + ".dat";
          fp = gmx_ffopen(ffh, "w");
          for (int k = 0; k < interm_cross_mat_density_[cross_index_[i][j]][ii][0].size(); k++)
          {
            fprintf(fp, "%lf", density_bins_[k]);
            for (int jj = 0; jj < natmol2_[j]; jj++)
            {
              fprintf(fp, " %lf", interm_cross_mat_density_[cross_index_[i][j]][ii][jj][k]);
            }
            fprintf(fp, "\n");
          }
          gmx_ffclose(fp);
        }
      }
    }
  }

  for (int i = 0; i < natmol2_.size(); i++)
  {
    FILE *fp = nullptr;
    std::string inter_file_name(outfile_inter_);
    std::size_t found = inter_file_name.find_last_of(".");
    fp = gmx_ffopen(inter_file_name.insert(found, "_" + std::to_string(i + 1) + "_" + std::to_string(i + 1)), "w");
    for (int ii = 0; ii < natmol2_[i]; ii++)
    {
      for (int jj = 0; jj < natmol2_[i]; jj++)
      {
        double dx = cutoff_ / static_cast<double>(interm_same_mat_density_[i][ii][jj].size());
        double dm = calc_mean(interm_same_mat_density_[i][ii][jj], dx);
        double prob = calc_prob(interm_same_mat_density_[i][ii][jj], dx);
        fprintf(fp, "%4i %4i %4i %4i %9.6lf %9.6lf\n", i + 1, ii + 1, i + 1, jj + 1, dm, prob);
      }
    }
    gmx_ffclose(fp);
    std::string intra_file_name(outfile_intra_);
    found = intra_file_name.find_last_of(".");
    fp = gmx_ffopen(intra_file_name.insert(found, "_" + std::to_string(i + 1) + "_" + std::to_string(i + 1)), "w");
    for (int ii = 0; ii < natmol2_[i]; ii++)
    {
      for (int jj = 0; jj < natmol2_[i]; jj++)
      {
        double dx = cutoff_ / static_cast<double>(intram_mat_density_[i][ii][jj].size());
        double dm = calc_mean(intram_mat_density_[i][ii][jj], dx);
        double prob = calc_prob(intram_mat_density_[i][ii][jj], dx);
        fprintf(fp, "%4i %4i %4i %4i %9.6lf %9.6lf\n", i + 1, ii + 1, i + 1, jj + 1, dm, prob);
      }
    }
    gmx_ffclose(fp);
    for (int j = i + 1; j < natmol2_.size(); j++)
    {
      std::string inter_c_file_name(outfile_inter_);
      found = inter_c_file_name.find_last_of(".");
      fp = gmx_ffopen(inter_c_file_name.insert(found, "_" + std::to_string(i + 1) + "_" + std::to_string(j + 1)), "w");
      for (int ii = 0; ii < natmol2_[i]; ii++)
      {
        for (int jj = 0; jj < natmol2_[j]; jj++)
        {
          double dx = cutoff_ / static_cast<double>(interm_cross_mat_density_[cross_index_[i][j]][ii][jj].size());
          double dm = calc_mean(interm_cross_mat_density_[cross_index_[i][j]][ii][jj], dx);
          double prob = calc_prob(interm_cross_mat_density_[cross_index_[i][j]][ii][jj], dx);
          fprintf(fp, "%4i %4i %4i %4i %9.6lf %9.6lf\n", i + 1, ii + 1, j + 1, jj + 1, dm, prob);
        }
      }
      gmx_ffclose(fp);
    }
  }
}

} // namespace

const char CMDataInfo::name[] = "cmdata";
const char CMDataInfo::shortDescription[] = "Calculate contact data";

TrajectoryAnalysisModulePointer CMDataInfo::create()
{
  return TrajectoryAnalysisModulePointer(new CMData);
}

} // namespace analysismodules

} // namespace gmx
