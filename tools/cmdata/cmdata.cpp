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

namespace gmx
{
 
namespace analysismodules
{
 
namespace
{

class Semaphore
{
private:
  std::mutex mut;
  std::condition_variable cv;
  std::uint64_t counter;

public:
  Semaphore( std::uint64_t counter = 0 ) : counter(counter), cv(std::condition_variable()), mut(std::mutex()) {}
  void acquire()
  {
    std::unique_lock<std::mutex> lock(mut);
    cv.wait(lock, [this](){ return counter > 0; });
    --counter;
    lock.unlock();
  }

  void release()
  {
    std::unique_lock<std::mutex> lock(mut);
    ++counter;
    lock.unlock();
    cv.notify_one();
  }

  void set_counter( std::uint64_t c )
  {
    counter = c;
  }

  std::uint64_t get_counter() const
  {
    return counter;
  }
};

using ftype_intra_ = void(
  int, std::size_t, std::size_t, double, int, const std::vector<int> &, const std::vector<int> &,
  const std::vector<double> &, const std::vector<double> &, std::vector<std::vector<std::mutex>> &, 
  std::vector<std::vector<std::vector<std::vector<double>>>> &
);
using ftype_same_ = void(
  int, std::size_t, std::size_t, std::size_t, double, const std::vector<int> &, const std::vector<int> &,
  const std::vector<double> &, std::vector<std::vector<std::mutex>> &, std::vector<std::vector<double>> &,
  std::vector<std::vector<std::vector<std::vector<double>>>> &
);
using ftype_cross_ = void(
  int, int, std::size_t, std::size_t, std::size_t, std::size_t, double, const std::vector<int> &,
  const std::vector<int> &, const std::vector<std::vector<int>> &, const std::vector<double> &,
  std::vector<std::vector<std::mutex>> &, std::vector<std::vector<double>> &,
  std::vector<std::vector<std::vector<std::vector<double>>>> &
);

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
  int nskip_;
  rvec *xcm_ = nullptr;
  gmx_mtop_t *mtop_;
  std::vector<int> mol_id_;
  std::string outfile_inter_;
  std::string outfile_intra_;
  std::vector<double> inv_num_mol_;
  std::vector<double> inv_num_mol_unique_;
  std::string sym_file_path_;
  std::vector<std::vector<std::vector<int>>> equivalence_list_;
  bool list_sym_;
  std::string list_sym_path_;
  std::vector<int> num_mol_unique_;

  std::vector<int> natmol2_;
  int nindex_;
  gmx::RangePartitioning mols_;
  std::vector<std::vector<int>> cross_index_;
  std::vector<t_atoms> molecules_;
  std::vector<double> density_bins_;
  std::size_t n_bins_;
  int num_threads_;

  double mcut2_;
  double cut_sig_2_;
  double dx_;
  std::vector<std::vector<std::vector<std::vector<double>>>> interm_same_mat_density_;
  std::vector<std::vector<std::vector<std::vector<double>>>> interm_cross_mat_density_;
  std::vector<std::vector<std::vector<std::vector<double>>>> intram_mat_density_;
  std::vector<std::vector<std::vector<std::vector<double>>>> interm_same_maxcdf_mol_;
  std::vector<std::vector<std::vector<std::vector<double>>>> interm_cross_maxcdf_mol_;

  // temporary containers for maxcdf operations
  std::vector<std::vector<double>> frame_same_mat_;
  std::vector<std::vector<double>> frame_cross_mat_;
  std::vector<std::vector<std::mutex>> frame_same_mutex_;
  std::vector<std::vector<std::mutex>> frame_cross_mutex_;
  std::vector<std::thread> threads_;
  std::vector<std::thread> mol_threads_;
  Semaphore semaphore_;

  // mode selection, booleans and functions
  std::string mode_;
  bool intra_ = false, same_ = false, cross_ = false;
  std::function<ftype_intra_> f_intra_mol_;
  std::function<ftype_same_> f_inter_mol_same_;
  std::function<ftype_cross_> f_inter_mol_cross_;
};

CMData::CMData() : sym_file_path_("") {}

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
  options->addOption(StringOption("ointra")
                    .store(&outfile_intra_)
                    .defaultValue("intramat.ndx")
                    .description("Output of the intra-molecular contacts"));
  options->addOption(StringOption("ointer")
                    .store(&outfile_inter_)
                    .defaultValue("intramat.ndx")
                    .description("Output of the intra-molecular contacts"));
  options->addOption(BooleanOption("histo")
                    .store(&histo_)
                    .defaultValue(false)
                    .description("Set to true to output histograms"));
  options->addOption(IntegerOption("nskip")
                    .store(&nskip_)
                    .defaultValue(0)
                    .description("Use every nskip-th frame"));
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
                    .defaultValue(false)
                    .description("Write symmetry list"));
  options->addOption(IntegerOption("num_threads")
                    .store(&num_threads_)
                    .defaultValue(1)
                    .description("Number of threads to run task on"));
  options->addOption(StringOption("mode")
                    .store(&mode_)
                    .defaultValue("intra+same+cross")
                    .description("Select the mode of analysis from intra, same or cross (combine with +)"));

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
  return prob;
}

static inline int n_bins(const double cut, const double factor = 4.0)
{
  return cut / (0.01 / factor);
}

static std::size_t offset_same( std::size_t i, int mol_i, std::size_t a_i, std::size_t a_j, const std::vector<int> &natmol2 )
{
  return static_cast<std::size_t>(mol_i) * (natmol2[i] * natmol2[i]) + a_i * natmol2[i] + a_j;
}

static std::size_t offset_same( std::size_t i, std::size_t mol_i, std::size_t a_i, std::size_t a_j, const std::vector<int> &natmol2 )
{
  return mol_i * (natmol2[i] * natmol2[i]) + a_i * natmol2[i] + a_j;
}

static std::size_t mutex_access( std::size_t i, std::size_t a_i, std::size_t a_j, const std::vector<int> &natmol2 )
{
  return a_i * natmol2[i] + a_j;
}

static std::size_t offset_cross( 
  std::size_t i, std::size_t j, std::size_t mol_i, std::size_t mol_j,
  std::size_t a_i, std::size_t a_j, const std::vector<int> &natmol2
)
{
  return mol_i * (mol_j * natmol2[i] * natmol2[j]) + mol_j * (natmol2[i] * natmol2[j]) + a_i * natmol2[j] + a_j; 
}

static void intra_mol_routine( 
  int i, std::size_t a_i, std::size_t a_j, double dx2, int nsym, const std::vector<int> &mol_id_,
  const std::vector<int> &natmol2_, const std::vector<double> &density_bins_,
  const std::vector<double> &inv_num_mol_, std::vector<std::vector<std::mutex>> &frame_same_mutex_, 
  std::vector<std::vector<std::vector<std::vector<double>>>> &intram_mat_density_
)
{
  std::size_t same_mutex_index = mutex_access(mol_id_[i], a_i, a_j, natmol2_);
  std::unique_lock lock(frame_same_mutex_[mol_id_[i]][same_mutex_index]);
  kernel_density_estimator(std::begin(intram_mat_density_[mol_id_[i]][a_i][a_j]), density_bins_, std::sqrt(dx2), inv_num_mol_[i]/nsym);
  lock.unlock();
}

static void inter_mol_same_routine(
  int i, std::size_t mol_i, std::size_t a_i, std::size_t a_j, double dx2, const std::vector<int> &mol_id_,
  const std::vector<int> &natmol2_, const std::vector<double> &density_bins_,
  std::vector<std::vector<std::mutex>> &frame_same_mutex_, std::vector<std::vector<double>> &frame_same_mat_,
  std::vector<std::vector<std::vector<std::vector<double>>>> &interm_same_mat_density_
)
{
  double dist = std::sqrt(dx2);
  std::size_t same_access_index = offset_same(mol_id_[i], mol_i, a_i, a_j, natmol2_);
  std::size_t same_mutex_index = mutex_access(mol_id_[i], a_i, a_j, natmol2_);
  std::unique_lock lock(frame_same_mutex_[mol_id_[i]][same_mutex_index]);
  kernel_density_estimator(std::begin(interm_same_mat_density_[mol_id_[i]][a_i][a_j]), density_bins_, dist, 1.0);
  frame_same_mat_[mol_id_[i]][same_access_index] = std::min(dist, frame_same_mat_[mol_id_[i]][same_access_index]);
  lock.unlock();
}

static void inter_mol_cross_routine(
  int i, int j, std::size_t mol_i, std::size_t mol_j, std::size_t a_i, std::size_t a_j, double dx2, const std::vector<int> &mol_id_,
  const std::vector<int> &natmol2_, const std::vector<std::vector<int>> &cross_index_, const std::vector<double> &density_bins_, 
  std::vector<std::vector<std::mutex>> &frame_cross_mutex_, std::vector<std::vector<double>> &frame_cross_mat_,
  std::vector<std::vector<std::vector<std::vector<double>>>> &interm_cross_mat_density_
)
{
  double dist = std::sqrt(dx2);
  std::size_t cross_access_index = offset_cross(mol_id_[i], mol_id_[j], mol_i, mol_j, a_i, a_j, natmol2_);
  std::size_t cross_mutex_index = mutex_access(mol_id_[j], a_i, a_j, natmol2_);
  std::unique_lock lock(frame_cross_mutex_[cross_index_[mol_id_[i]][mol_id_[j]]][cross_mutex_index]);
  kernel_density_estimator(std::begin(interm_cross_mat_density_[cross_index_[mol_id_[i]][mol_id_[j]]][a_i][a_j]), density_bins_, dist, 1.0);
  frame_cross_mat_[cross_index_[mol_id_[i]][mol_id_[j]]][cross_access_index] = std::min(dist, frame_cross_mat_[cross_index_[mol_id_[i]][mol_id_[j]]][cross_access_index]);
  lock.unlock();
}

void CMData::initAnalysis(const TrajectoryAnalysisSettings &settings, const TopologyInformation &top)
{
  /* Checking inputs */
  GMX_RELEASE_ASSERT(num_threads_ > 0, "num_threads cannot be less than 1\n");
  if ( num_threads_ > std::thread::hardware_concurrency() )
  {
    num_threads_ = std::thread::hardware_concurrency();
    printf("Maximum thread number surpassed. Scaling num_threads down accordingly\n");
  }
  printf("Running operations on %i threads\n", num_threads_);
  threads_.resize(num_threads_);

  GMX_RELEASE_ASSERT(nskip_ >= 0, "nskip needs to be greater than or equal to 0\n");

  n_x_ = 0;
  nframe_ = 0;
  mtop_ = top.mtop();
  mols_ = gmx_mtop_molecules(*top.mtop());
  // number of molecules
  nindex_ = mols_.numBlocks();

  mol_threads_.resize(nindex_);
  semaphore_.set_counter(std::min(num_threads_, nindex_));

  // set up mode selection
  f_intra_mol_ = [](int, std::size_t, std::size_t, double, int, const std::vector<int> &,
  const std::vector<int> &, const std::vector<double> &, const std::vector<double> &, 
  std::vector<std::vector<std::mutex>> &, std::vector<std::vector<std::vector<std::vector<double>>>> &){};
  f_inter_mol_same_ = [](int, std::size_t, std::size_t, std::size_t, double, const std::vector<int> &, const std::vector<int> &,
  const std::vector<double> &, std::vector<std::vector<std::mutex>> &, std::vector<std::vector<double>> &,
  std::vector<std::vector<std::vector<std::vector<double>>>> &){};
  f_inter_mol_cross_ = [](int, int, std::size_t, std::size_t, std::size_t, std::size_t, double, const std::vector<int> &,
  const std::vector<int> &, const std::vector<std::vector<int>> &, const std::vector<double> &, std::vector<std::vector<std::mutex>> &,
  const std::vector<std::vector<double>> &, std::vector<std::vector<std::vector<std::vector<double>>>> &){};

  printf("Evaluating mode selection:\n");
  std::string tmp_mode;
  std::stringstream modestream{ mode_ };
  while (std::getline(modestream, tmp_mode, '+'))
  {
    printf(" - found %s", tmp_mode.c_str());
    if ( tmp_mode == std::string("intra") )
    {
      intra_ = true;
      f_intra_mol_ = intra_mol_routine;
      printf(" :: activating intramat calculations\n");
    }
    else if ( tmp_mode == std::string("same") )
    {
      same_ = true;
      f_inter_mol_same_ = inter_mol_same_routine;
      printf(" :: activating intermat same calculations\n");
    }
    else if ( tmp_mode == std::string("cross") )
    {
      cross_ = true;
      f_inter_mol_cross_ = inter_mol_cross_routine;
      printf(" :: activating intermat cross calculations\n");
    }
    else printf(" :: ignoring keyword\n");
  }

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
  num_mol_unique_.push_back(num_mol[num_unique_molecules]);
  inv_num_mol_unique_.push_back(1. / (static_cast<double>(num_mol[num_unique_molecules])));
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
      inv_num_mol_unique_.push_back(1. / (static_cast<double>(num_mol[num_unique_molecules])));
      num_mol_unique_.push_back(num_mol[num_unique_molecules]);
    }
    mol_id_.push_back(num_unique_molecules);
    inv_num_mol_.push_back(1. / static_cast<double>(num_mol[num_unique_molecules]));
  }

  printf("number of different molecules %lu\n", natmol2_.size());
  for(std::size_t i=0; i<natmol2_.size();i++) printf("mol %lu num %u size %u\n", i, num_mol[i], natmol2_[i]);

  if (same_)
  {
    interm_same_mat_density_.resize(natmol2_.size());
    interm_same_maxcdf_mol_.resize(natmol2_.size());
  }
  if (cross_)
  {
    interm_cross_mat_density_.resize((natmol2_.size() * (natmol2_.size() - 1)) / 2);
    interm_cross_maxcdf_mol_.resize((natmol2_.size() * (natmol2_.size() - 1)) / 2);
  }
  if (intra_) intram_mat_density_.resize(natmol2_.size());

  density_bins_.resize(n_bins(cutoff_));
  for (std::size_t i = 0; i < density_bins_.size(); i++)
    density_bins_[i] = cutoff_ / static_cast<double>(density_bins_.size()) * static_cast<double>(i) + cutoff_ / static_cast<double>(density_bins_.size() * 2);

  int cross_count = 0;
  if (cross_) cross_index_.resize(natmol2_.size(), std::vector<int>(natmol2_.size(), 0));
  for ( std::size_t i = 0; i < natmol2_.size(); i++ )
  {
    if (same_)
    {
      interm_same_mat_density_[i].resize(natmol2_[i], std::vector<std::vector<double>>(natmol2_[i], std::vector<double>(n_bins(cutoff_), 0)));
      interm_same_maxcdf_mol_[i].resize(natmol2_[i], std::vector<std::vector<double>>(natmol2_[i], std::vector<double>(n_bins(cutoff_), 0)));
    }
    if (intra_) intram_mat_density_[i].resize(natmol2_[i], std::vector<std::vector<double>>(natmol2_[i], std::vector<double>(n_bins(cutoff_), 0)));
    for ( std::size_t j = i + 1; j < natmol2_.size() && cross_; j++ )
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
      for ( int j = 0; j < equivalence_list_[i].size(); j++ )
      {
        sym_list_file << "atom " << j << ":";
        for ( int k = 0; k < equivalence_list_[i][j].size(); k++ )
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
  dx_ = cutoff_ / static_cast<double>(n_bins_);

  mcut2_ = mol_cutoff_ * mol_cutoff_;
  cut_sig_2_ = (cutoff_ + 0.02) * (cutoff_ + 0.02);
  snew(xcm_, nindex_);

  if (same_) frame_same_mat_.resize(natmol2_.size());
  if (intra_ | same_) frame_same_mutex_.resize(natmol2_.size());
  if (cross_) frame_cross_mat_.resize(cross_index_.size());
  if (cross_) frame_cross_mutex_.resize(cross_index_.size());
  std::size_t sum_cross_mol_sizes = 0;
  for ( std::size_t i = 0; i < natmol2_.size(); i++ )
  {
    if (same_) frame_same_mat_[i].resize(natmol2_[i] *  natmol2_[i] * num_mol_unique_[i], 0);
    if (intra_ | same_) frame_same_mutex_[i] = std::vector<std::mutex>(natmol2_[i] *  natmol2_[i]);
    for ( std::size_t j = i+1; j < natmol2_.size() && cross_; j++ )
    {
      frame_cross_mat_[cross_index_[i][j]].resize(natmol2_[i] * natmol2_[j] * num_mol_unique_[i] * num_mol_unique_[j], 0);
      frame_cross_mutex_[cross_index_[i][j]] = std::vector<std::mutex>(natmol2_[i] * natmol2_[j]);
    }
  }

  printf("Finished preprocessing. Starting frame-by-frame analysis.\n");
}

static void mindist_same(
  std::size_t start_mti_same, std::size_t start_im_same, std::size_t start_i_same, 
  std::size_t start_j_same, long int n_loop_operations_same, const std::vector<double> &density_bins,
  const std::vector<int> &num_mol_unique, const std::vector<int> &natmol2,
  const std::vector<std::vector<double>> &frame_same_mat,
  std::vector<std::vector<std::mutex>> &frame_same_mutex,
  std::vector<std::vector<std::vector<std::vector<double>>>> &interm_same_maxcdf_mol
)
{
  bool first_im_same = true, first_i_same = true, first_j_same = true;
  unsigned counter_same = 0;
  for (std::size_t mt_i = start_mti_same; mt_i < natmol2.size(); mt_i++)
  {
    for (std::size_t im = (first_im_same) ? start_im_same : 0; im < num_mol_unique[mt_i]; im++)
    {
      for (std::size_t i = (first_i_same) ? start_i_same : 0; i < natmol2[mt_i]; i++)
      {
        for (std::size_t j = (first_j_same) ? start_j_same : i; j < natmol2[mt_i]; j++)
        {
          std::size_t offset = offset_same(mt_i, im, i, j, natmol2);
          std::size_t mutex_j = mutex_access(mt_i, i, j, natmol2);
          double mindist = frame_same_mat[mt_i][offset];

          std::unique_lock<std::mutex> lock(frame_same_mutex[mt_i][mutex_j]);
          kernel_density_estimator(std::begin(interm_same_maxcdf_mol[mt_i][i][j]), density_bins, mindist, 1.0);
          interm_same_maxcdf_mol[mt_i][j][i] = interm_same_maxcdf_mol[mt_i][i][j];
          lock.unlock();

          ++counter_same;
          if ( counter_same == n_loop_operations_same ) return;
        }
        first_j_same = false;
      }
      first_i_same = false;
    }
    first_im_same = false;
  }
}

static void mindist_cross(
  std::size_t start_mti_cross, std::size_t start_mtj_cross, std::size_t start_im_cross, std::size_t start_jm_cross, 
  std::size_t start_i_cross, std::size_t start_j_cross, int n_loop_operations_cross, 
  const std::vector<int> &natmol2, const std::vector<std::vector<int>> &cross_index,
  const std::vector<double> &density_bins, const std::vector<int> &num_mol_unique,
  const std::vector<std::vector<double>> &frame_cross_mat,
  std::vector<std::vector<std::mutex>> &frame_cross_mutex,
  std::vector<std::vector<std::vector<std::vector<double>>>> &interm_cross_maxcdf_mol
)
{
  bool first_im_cross = true, first_mtj_cross = true, first_jm_cross = true, first_i_cross = true, first_j_cross = true;
  unsigned counter_cross = 0;
  
  for (std::size_t mt_i = start_mti_cross; mt_i < natmol2.size(); mt_i++)
  {
    for (std::size_t mt_j = (first_mtj_cross) ? start_mtj_cross : mt_i + 1; mt_j < natmol2.size(); mt_j++)
    {
      for (std::size_t im = (first_im_cross) ? start_im_cross : 0; im < num_mol_unique[mt_i]; im++)
      {
        for (std::size_t jm = (first_jm_cross) ? start_jm_cross : 0; jm < num_mol_unique[mt_j]; jm++)
        {
          for (std::size_t i = (first_i_cross) ? start_i_cross : 0; i < natmol2[mt_i]; i++)
          {
            for (std::size_t j = (first_j_cross) ? start_j_cross : 0; j < natmol2[mt_j]; j++)
            {
              std::size_t offset = offset_cross(mt_i, mt_j, im, jm, i, j, natmol2);
              std::size_t mutex_j = mutex_access(mt_j, i, j, natmol2);
              double mindist = frame_cross_mat[cross_index[mt_i][mt_j]][offset];
              
              std::unique_lock<std::mutex> lock(frame_cross_mutex[cross_index[mt_i][mt_j]][mutex_j]);
              kernel_density_estimator(std::begin(interm_cross_maxcdf_mol[cross_index[mt_i][mt_j]][i][j]), density_bins, mindist, 1.0);
              lock.unlock();

              ++counter_cross;
              if ( counter_cross == n_loop_operations_cross ) return;
            }
            first_j_cross = false;
          }
          first_i_cross = false;
        }
        first_jm_cross = false;
      }
      first_mtj_cross = false;
    }
    first_im_cross = false;
  }
}

static void mindist_kernel(
  const std::vector<int> &natmol2,          // common parameters
  const std::vector<double> &density_bins,
  const std::vector<int> &num_mol_unique,
  std::size_t start_mti_same,               // same parameters
  std::size_t start_im_same,
  std::size_t start_i_same,
  std::size_t start_j_same,
  long int n_loop_operations_same,
  const std::vector<std::vector<double>> &frame_same_mat,
  std::vector<std::vector<std::mutex>> &frame_same_mutex, 
  std::vector<std::vector<std::vector<std::vector<double>>>> &interm_same_maxcdf_mol,
  std::size_t start_mti_cross,              // cross parameters
  std::size_t start_mtj_cross,
  std::size_t start_im_cross,
  std::size_t start_jm_cross,
  std::size_t start_i_cross,
  std::size_t start_j_cross,
  int n_loop_operations_cross,
  const std::vector<std::vector<int>> &cross_index, 
  const std::vector<std::vector<double>> &frame_cross_mat,
  std::vector<std::vector<std::mutex>> &frame_cross_mutex,
  std::vector<std::vector<std::vector<std::vector<double>>>> &interm_cross_maxcdf_mol
)
{
  if ( n_loop_operations_same != 0 )
  {
    mindist_same(
      start_mti_same, start_im_same, start_i_same, start_j_same, n_loop_operations_same, density_bins,
      num_mol_unique, natmol2, frame_same_mat, frame_same_mutex, interm_same_maxcdf_mol
    );
  }

  if ( n_loop_operations_cross != 0 )
  {
    mindist_cross(
      start_mti_cross, start_mtj_cross, start_im_cross, start_jm_cross, start_i_cross, start_j_cross, n_loop_operations_cross, 
      natmol2, cross_index, density_bins, num_mol_unique, frame_cross_mat, frame_cross_mutex, interm_cross_maxcdf_mol
    );
  }
}

static void molecule_routine(
  const int i, const int nindex_, t_pbc *pbc, rvec *x, const std::vector<double> &inv_num_mol_, const double cut_sig_2_, 
  const std::vector<int> &natmol2_, const std::vector<int> &num_mol_unique_, const std::vector<int> &mol_id_, 
  const std::vector<std::vector<int>> &cross_index_, const std::vector<double> &density_bins_, const double mcut2_, 
  rvec *xcm_, const gmx::RangePartitioning &mols_, gmx_mtop_t *mtop_, const std::vector<std::vector<std::vector<int>>> &equivalence_list_,
  std::vector<std::vector<double>> &frame_same_mat_, std::vector<std::vector<std::mutex>> &frame_same_mutex_,
  std::vector<std::vector<std::vector<std::vector<double>>>> &intram_mat_density_,
  std::vector<std::vector<std::vector<std::vector<double>>>> &interm_same_mat_density_,
  std::vector<std::vector<double>> &frame_cross_mat_, std::vector<std::vector<std::mutex>> &frame_cross_mutex_,
  std::vector<std::vector<std::vector<std::vector<double>>>> &interm_cross_mat_density_, Semaphore &semaphore_,
  const std::function<ftype_intra_> &f_intra_mol_, const std::function<ftype_same_> &f_inter_mol_same_, const std::function<ftype_cross_> &f_inter_mol_cross_
)
{
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
  semaphore_.acquire();
  /* Loop over molecules  */
  for (int j = 0; j < nindex_; j++)
  {
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
    GMX_RELEASE_ASSERT(mols_.numBlocks() > 0, "Cannot access index[] from empty mols");
    /* cycle over the atoms of a molecule i */
    for (std::size_t ii = mols_.block(i).begin(); ii < mols_.block(i).end(); ii++)
    {
      std::size_t a_j = 0;
      mtopGetAtomAndResidueName(*mtop_, ii, &molb, &atomname, nullptr, nullptr, nullptr);
      if (atomname[0] == 'H')
      {
        a_i++;
        continue;
      }
      /* cycle over the atoms of a molecule j */
      for (std::size_t jj = mols_.block(j).begin(); jj < mols_.block(j).end(); jj++)
      {
        mtopGetAtomAndResidueName(*mtop_, jj, &molb, &atomname, nullptr, nullptr, nullptr);
        if (atomname[0] == 'H')
        {
          a_j++;
          continue;
        }
        // check for chemical equivalence
        for (std::size_t eq_i = 0; eq_i < equivalence_list_[mol_id_[i]][a_i].size(); eq_i++)
        {
          for (std::size_t eq_j = 0; eq_j < equivalence_list_[mol_id_[j]][a_j].size(); eq_j++)
          {
            // get molecule-wise atom index considering equivalence
            std::size_t eqa_i  = equivalence_list_[mol_id_[i]][a_i][eq_i];             // molecule-wise equivalence index i
            std::size_t geqa_i = ii + (eqa_i - equivalence_list_[mol_id_[i]][a_i][0]); // global equivalence index i
            std::size_t eqa_j  = equivalence_list_[mol_id_[j]][a_j][eq_j];             // molecule-wise equivalence index j
            std::size_t geqa_j = jj + (eqa_j - equivalence_list_[mol_id_[j]][a_j][0]); // global equivalence index j
            std::size_t delta  = eqa_i - eqa_j;
            double nsym = static_cast<double>(equivalence_list_[mol_id_[i]][a_i].size()*equivalence_list_[mol_id_[i]][a_j].size());
            rvec sym_dx;
            if (pbc != nullptr) pbc_dx(pbc, x[geqa_i], x[geqa_j], sym_dx);
            else rvec_sub(x[geqa_i], x[geqa_j], sym_dx);
            double dx2 = iprod(sym_dx, sym_dx);
            if(i==j) 
            {
              if (dx2 < cut_sig_2_)
              {
                f_intra_mol_(i, a_i, a_j, dx2, nsym, mol_id_, natmol2_, density_bins_, inv_num_mol_, frame_same_mutex_, intram_mat_density_);
              }
            }
            else
            {
              if(mol_id_[i]==mol_id_[j])
              { // inter same molecule specie
                if (dx2 < cut_sig_2_)
                {
                  f_inter_mol_same_(
                    i, mol_i, a_i, a_j, dx2, mol_id_, natmol2_, density_bins_, frame_same_mutex_, frame_same_mat_, interm_same_mat_density_
                  );
                }
                if(delta!=0.) {
                  // this is to account for inversion atom/molecule
                  if (pbc != nullptr) pbc_dx(pbc, x[geqa_i-delta], x[geqa_j+delta], sym_dx);
                  else rvec_sub(x[geqa_i-delta], x[geqa_j+delta], sym_dx);
                  dx2 = iprod(sym_dx, sym_dx);
                  if (dx2 < cut_sig_2_)
                  {
                    f_inter_mol_same_(
                      i, mol_i, a_i, a_j, dx2, mol_id_, natmol2_, density_bins_, frame_same_mutex_, frame_same_mat_, interm_same_mat_density_
                    );
                  }
                }
              } 
              else
              { // inter cross molecule species
                if (dx2 < cut_sig_2_)
                {
                  f_inter_mol_cross_(
                    i, j, mol_i, mol_j, a_i, a_j, dx2, mol_id_, natmol2_, cross_index_, density_bins_, frame_cross_mutex_, frame_cross_mat_, interm_cross_mat_density_
                  );
                }
              }
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

void CMData::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc, TrajectoryAnalysisModuleData *pdata)
{
  if ((nskip_ == 0) || ((nskip_ > 0) && ((frnr % nskip_) == 0)))
  {
    rvec *x = fr.x;
    /* resetting the per frame vector to zero */
    for ( std::size_t i = 0; i < frame_same_mat_.size(); i++ )
    {
      #pragma omp parallel for num_threads(std::min(num_threads_, static_cast<int>(frame_same_mat_[i].size())))
      for ( std::size_t j = 0; j < frame_same_mat_[i].size(); j++ ) frame_same_mat_[i][j] = 100.;
    }
    for ( std::size_t i = 0; i < frame_cross_mat_.size(); i++ )
    {
      #pragma omp parallel for num_threads(std::min(num_threads_, static_cast<int>(frame_cross_mat_[i].size())))
      for ( std::size_t j = 0; j < frame_cross_mat_[i].size(); j++ ) frame_cross_mat_[i][j] = 100.;
    }

    /* start loop for each molecule */
    for (int i = 0; i < nindex_; i++)
    {
      /* start molecule thread*/
      mol_threads_[i] = std::thread(molecule_routine, i, nindex_, pbc, x, std::cref(inv_num_mol_), 
      cut_sig_2_, std::cref(natmol2_), std::cref(num_mol_unique_), std::cref(mol_id_), std::cref(cross_index_),
      std::cref(density_bins_), mcut2_, xcm_, mols_, mtop_, std::cref(equivalence_list_), std::ref(frame_same_mat_),
      std::ref(frame_same_mutex_), std::ref(intram_mat_density_), std::ref(interm_same_mat_density_), std::ref(frame_cross_mat_), 
      std::ref(frame_cross_mutex_), std::ref(interm_cross_mat_density_), std::ref(semaphore_), std::cref(f_intra_mol_),
      std::cref(f_inter_mol_same_), std::cref(f_inter_mol_cross_));
      /* end molecule thread */
    }
    /* join molecule threads */
    for ( auto &thread : mol_threads_ ) thread.join();

    /* calculate the mindist accumulation indices */
    std::size_t num_ops_same = 0;
    for ( std::size_t im = 0; im < natmol2_.size(); im++ ) num_ops_same += num_mol_unique_[im] * ( natmol2_[im] * ( natmol2_[im] + 1 ) ) / 2;
    int n_per_thread_same = (same_) ? num_ops_same / num_threads_  : 0;
    int n_threads_same_uneven = (same_) ? num_ops_same % num_threads_ : 0;
    std::size_t start_mti_same = 0, start_im_same = 0, end_mti_same = 1, end_im_same = 1; 
    std::size_t start_i_same = 0, start_j_same = 0, end_i_same = 0, end_j_same = 0;
    int num_ops_cross = 0;
    for ( std::size_t im = 0; im < natmol2_.size(); im++ )
    {
      for ( std::size_t jm = im + 1; jm < natmol2_.size(); jm++ )
      {
        num_ops_cross += num_mol_unique_[im] * natmol2_[im] * num_mol_unique_[jm] * natmol2_[jm];
      }
    }
    int n_per_thread_cross = (cross_) ? num_ops_cross / num_threads_ : 0;
    int n_threads_cross_uneven = (cross_) ? num_ops_cross % num_threads_ : 0;

    std::size_t start_mti_cross = 0, start_mtj_cross = 1, start_im_cross = 0, start_jm_cross = 0, start_i_cross = 0, start_j_cross = 0;
    std::size_t end_mti_cross = 1, end_mtj_cross = 2, end_im_cross = 1, end_jm_cross = 1, end_i_cross = 0, end_j_cross = 0;
    
    for ( int tid = 0; tid < num_threads_; tid++ )
    {
      /* calculate same indices */
      int n_loop_operations_same = n_per_thread_same + (tid < n_threads_same_uneven ? 1 : 0);
      long int n_loop_operations_total_same = n_loop_operations_same;
      while ( natmol2_[end_mti_same - 1] - static_cast<int>(end_j_same) <= n_loop_operations_same )
      {
        int sub_same = natmol2_[end_mti_same - 1] - static_cast<int>(end_j_same);
        n_loop_operations_same -= sub_same;
        end_i_same++;
        end_j_same = end_i_same;
        if ( static_cast<int>(end_j_same) == natmol2_[end_mti_same - 1] )
        {
          end_im_same++;
          end_i_same = 0;
          end_j_same = 0;
        }
        if ( static_cast<int>(end_im_same) == num_mol_unique_[end_mti_same - 1] )
        {
          end_mti_same++;
          end_im_same = 1;
          end_i_same = 0;
          end_j_same = 0;
        }
        if (n_loop_operations_same == 0) break;
      }
      end_j_same += n_loop_operations_same;  

      /* calculate cross indices */
      int n_loop_operations_total_cross = n_per_thread_cross + ( tid < n_threads_cross_uneven ? 1 : 0 );
      if ( natmol2_.size() > 1 && cross_ )
      {
        int n_loop_operations_cross = n_loop_operations_total_cross;
        while ( natmol2_[end_mti_cross-1] * natmol2_[end_mtj_cross-1] - (natmol2_[end_mtj_cross-1] * static_cast<int>(end_i_cross) + static_cast<int>(end_j_cross)) <= n_loop_operations_cross )
        {
          int sub_cross = natmol2_[end_mti_cross-1] * natmol2_[end_mtj_cross-1] - (natmol2_[end_mtj_cross-1] * static_cast<int>(end_i_cross) + static_cast<int>(end_j_cross));
          n_loop_operations_cross -= sub_cross;

          end_jm_cross++;
          end_i_cross = 0;
          end_j_cross = 0;

          // case jm is above max
          if ( end_jm_cross > num_mol_unique_[end_mtj_cross - 1] )
          {
            end_im_cross++;
            end_jm_cross = 1;
            end_i_cross = 0;
            end_j_cross = 0;
          }
          // case im is above max
          if ( end_im_cross > num_mol_unique_[end_mti_cross - 1] )
          {
            end_mtj_cross++;
            end_im_cross = 1;
            end_jm_cross = 1;
            end_i_cross = 0;
            end_j_cross = 0;
          }
          // case mtj is above max
          if ( end_mtj_cross > natmol2_[end_mtj_cross - 1] )
          {
            end_mti_cross++;
            end_mtj_cross = end_mti_cross + 1;
            end_im_cross = 1;
            end_jm_cross = 1;
            end_i_cross = 0;
            end_j_cross = 0;
          }
          if ( n_loop_operations_cross == 0 ) break;
          if ( end_mti_cross == natmol2_.size() - 1 ) break;
        }

        // calculate overhangs and add them
        end_i_cross += n_loop_operations_cross / natmol2_[end_mtj_cross-1];
        end_j_cross += n_loop_operations_cross % natmol2_[end_mtj_cross-1];
        end_i_cross += end_j_cross / natmol2_[end_mtj_cross-1];
        end_j_cross %= natmol2_[end_mtj_cross-1];
      }

      /* start thread */
      threads_[tid] = std::thread(
        mindist_kernel, std::cref(natmol2_), std::cref(density_bins_), std::cref(num_mol_unique_), 
        start_mti_same, start_im_same, start_i_same, start_j_same, n_loop_operations_total_same,
        std::cref(frame_same_mat_), std::ref(frame_same_mutex_), std::ref(interm_same_maxcdf_mol_),
        start_mti_cross, start_mtj_cross, start_im_cross, start_jm_cross, start_i_cross, start_j_cross,
        n_loop_operations_total_cross, std::cref(cross_index_), std::cref(frame_cross_mat_), 
        std::ref(frame_cross_mutex_), std::ref(interm_cross_maxcdf_mol_)
      );
      /* end thread */

      /* set new starts */
      start_mti_same = end_mti_same - 1;
      start_im_same = end_im_same - 1;
      start_i_same = end_i_same;
      start_j_same = end_j_same;

      start_mti_cross = end_mti_cross - 1;
      start_mtj_cross = end_mtj_cross - 1;
      start_im_cross = end_im_cross - 1;
      start_jm_cross = end_jm_cross - 1;
      start_i_cross = end_i_cross;
      start_j_cross = end_j_cross;
    }
    for ( auto &thread : threads_ ) thread.join();
    ++n_x_;
  }
  ++nframe_;
}

static void normalize_histo(
  std::size_t i, int ii, int jj, double norm, double inv_num_mol_same,
  std::vector<std::vector<std::vector<std::vector<double>>>> &data
)
{
  std::transform(std::begin(data[i][ii][jj]), std::end(data[i][ii][jj]), std::begin(data[i][ii][jj]),
                 [norm, inv_num_mol_same]( auto &c ){ return c * norm * inv_num_mol_same; }
                );
}

void CMData::finishAnalysis(int /*nframes*/)
{
  // normalisations
  double norm = 1. / n_x_;

  using ftype_norm = void(std::size_t, int, int, double, double, std::vector<std::vector<std::vector<std::vector<double>>>> &);
  auto f_empty = [](std::size_t, int, int, double, double, std::vector<std::vector<std::vector<std::vector<double>>>> &){};

  std::function<ftype_norm> normalize_intra = (intra_) ? normalize_histo : f_empty;
  std::function<ftype_norm> normalize_inter_same = (same_) ? normalize_histo : f_empty;
  std::function<ftype_norm> normalize_inter_cross = (cross_) ? normalize_histo : f_empty;

  for (std::size_t i = 0; i < natmol2_.size(); i++)
  {
    for (int ii = 0; ii < natmol2_[i]; ii++)
    {
      for (int jj = ii; jj < natmol2_[i]; jj++)
      {
        double inv_num_mol_same = inv_num_mol_unique_[i];
        normalize_inter_same(i, ii, jj, norm, inv_num_mol_same, interm_same_maxcdf_mol_);
        normalize_inter_same(i, ii, jj, norm, 1.0, interm_same_mat_density_);
        normalize_intra(i, ii, jj, norm, 1.0, intram_mat_density_);
        
        double sum = 0.0;
        for ( std::size_t k = (same_) ? 0 : n_bins(cutoff_); k < n_bins(cutoff_); k++ )
        {
          sum+= dx_ * interm_same_maxcdf_mol_[i][ii][jj][k];
          if ( sum > 1.0 ) sum=1.0;
          interm_same_maxcdf_mol_[i][ii][jj][k] = sum;
        }
        if (same_) interm_same_mat_density_[i][jj][ii] = interm_same_mat_density_[i][ii][jj];
        if (same_) interm_same_maxcdf_mol_[i][jj][ii] = interm_same_maxcdf_mol_[i][ii][jj];
        if (intra_) intram_mat_density_[i][jj][ii] = intram_mat_density_[i][ii][jj];
      }
    }
    for (std::size_t j = i + 1; j < natmol2_.size() & cross_; j++)
    {
      for (int ii = 0; ii < natmol2_[i]; ii++)
      {
        for (int jj = 0; jj < natmol2_[j]; jj++)
        {
          double inv_num_mol_cross = inv_num_mol_unique_[i];
          normalize_inter_cross(cross_index_[i][j], ii, jj, norm, 1.0, interm_cross_mat_density_);
          normalize_inter_cross(cross_index_[i][j], ii, jj, norm, inv_num_mol_cross, interm_cross_maxcdf_mol_);

          double sum = 0.0;
          for ( std::size_t k = (cross_) ? 0 : n_bins(cutoff_); k < n_bins(cutoff_); k++ )
          {
            sum += dx_ * interm_cross_maxcdf_mol_[cross_index_[i][j]][ii][jj][k];
            if ( sum > 1.0 ) sum = 1.0;
            interm_cross_maxcdf_mol_[cross_index_[i][j]][ii][jj][k] = sum;
          }
        }
      }
    }
  }
}

static void f_write_intra( 
  std::size_t i, int ii, const std::vector<double> &density_bins, const std::vector<int> &natmol2,
  const std::vector<std::vector<std::vector<std::vector<double>>>> &intram_mat_density
)
{
  FILE *fp_intra = nullptr;

  std::string ffh_intra = "intra_mol_" + std::to_string(i + 1) + "_" + std::to_string(i + 1) + "_aa_" + std::to_string(ii + 1) + ".dat";
  fp_intra = gmx_ffopen(ffh_intra, "w");
  for (int k = 0; k < density_bins.size(); k++)
  {
    fprintf(fp_intra, "%lf", density_bins[k]);
    for (int jj = 0; jj < natmol2[i]; jj++)
    {
      fprintf(fp_intra, " %lf", intram_mat_density[i][ii][jj][k]);
    }
    fprintf(fp_intra, "\n");
  }

  gmx_ffclose(fp_intra);
}

static void f_write_inter_same(
  std::size_t i, int ii, const std::vector<double> &density_bins, const std::vector<int> &natmol2,
  const std::vector<std::vector<std::vector<std::vector<double>>>> &interm_same_mat_density,
  const std::vector<std::vector<std::vector<std::vector<double>>>> &interm_same_maxcdf_mol
)
{
  FILE *fp_inter = nullptr;
  FILE *fp_inter_cum = nullptr;
  std::string ffh_inter = "inter_mol_" + std::to_string(i + 1) + "_" + std::to_string(i + 1) + "_aa_" + std::to_string(ii + 1) + ".dat";
  fp_inter = gmx_ffopen(ffh_inter, "w");
  std::string ffh_inter_cum = "inter_mol_c_" + std::to_string(i + 1) + "_" + std::to_string(i + 1) + "_aa_" + std::to_string(ii + 1) + ".dat";
  fp_inter_cum = gmx_ffopen(ffh_inter_cum, "w");
  for (int k = 0; k < density_bins.size(); k++)
  {
    fprintf(fp_inter, "%lf", density_bins[k]);
    fprintf(fp_inter_cum, "%lf", density_bins[k]);
    for (int jj = 0; jj < natmol2[i]; jj++)
    {
      fprintf(fp_inter, " %lf", interm_same_mat_density[i][ii][jj][k]);
      fprintf(fp_inter_cum, " %lf", interm_same_maxcdf_mol[i][ii][jj][k]);
    }
    fprintf(fp_inter, "\n");
    fprintf(fp_inter_cum, "\n");
  }
  gmx_ffclose(fp_inter);
  gmx_ffclose(fp_inter_cum);
}

static void f_write_inter_cross(
  std::size_t i, std::size_t j, int ii, const std::vector<double> &density_bins, const std::vector<int> &natmol2,
  const std::vector<std::vector<int>> &cross_index,
  const std::vector<std::vector<std::vector<std::vector<double>>>> &interm_cross_mat_density,
  const std::vector<std::vector<std::vector<std::vector<double>>>> &interm_cross_maxcdf_mol
)
{
  FILE *fp = nullptr;
  FILE *fp_cum = nullptr;
  std::string ffh = "inter_mol_" + std::to_string(i + 1) + "_" + std::to_string(j + 1) + "_aa_" + std::to_string(ii + 1) + ".dat";
  std::string ffh_cum = "inter_mol_c_" + std::to_string(i + 1) + "_" + std::to_string(j + 1) + "_aa_" + std::to_string(ii + 1) + ".dat";
  fp = gmx_ffopen(ffh, "w");
  fp_cum = gmx_ffopen(ffh_cum, "w");
  for (int k = 0; k < interm_cross_mat_density[cross_index[i][j]][ii][0].size(); k++)
  {
    fprintf(fp, "%lf", density_bins[k]);
    fprintf(fp_cum, "%lf", density_bins[k]);
    for (int jj = 0; jj < natmol2[j]; jj++)
    {
      fprintf(fp, " %lf", interm_cross_mat_density[cross_index[i][j]][ii][jj][k]);
      fprintf(fp_cum, " %lf", interm_cross_maxcdf_mol[cross_index[i][j]][ii][jj][k]);
    }
    fprintf(fp, "\n");
    fprintf(fp_cum, "\n");
  }
  gmx_ffclose(fp);
  gmx_ffclose(fp_cum);
}

void CMData::writeOutput()
{
  if (histo_)
  {
    printf("Writing data... ");
    std::function<void(
      std::size_t, int, const std::vector<double> &, const std::vector<int> &,
      const std::vector<std::vector<std::vector<std::vector<double>>>> &
    )> write_intra = [](
      std::size_t, int, const std::vector<double> &, const std::vector<int> &,
      const std::vector<std::vector<std::vector<std::vector<double>>>> &
    ){};

    std::function<void(
      std::size_t, int, const std::vector<double> &, const std::vector<int> &,
      const std::vector<std::vector<std::vector<std::vector<double>>>> &,
      const std::vector<std::vector<std::vector<std::vector<double>>>> &
    )> write_inter_same = [](
      std::size_t, int, const std::vector<double> &, const std::vector<int> &,
      const std::vector<std::vector<std::vector<std::vector<double>>>> &,
      const std::vector<std::vector<std::vector<std::vector<double>>>> &
    ){};

    std::function<void(
      std::size_t, std::size_t, int, const std::vector<double> &, const std::vector<int> &,
      const std::vector<std::vector<int>> &,
      const std::vector<std::vector<std::vector<std::vector<double>>>> &,
      const std::vector<std::vector<std::vector<std::vector<double>>>> &
    )> write_inter_cross = [](
      std::size_t, std::size_t, int, const std::vector<double> &, const std::vector<int> &,
      const std::vector<std::vector<int>> &,
      const std::vector<std::vector<std::vector<std::vector<double>>>> &,
      const std::vector<std::vector<std::vector<std::vector<double>>>> &
    ){};

    if (intra_) write_intra = f_write_intra;
    if (same_) write_inter_same = f_write_inter_same;
    if (cross_) write_inter_cross = f_write_inter_cross;
    for (std::size_t i = 0; i < natmol2_.size(); i++)
    {
      for (int ii = 0; ii < natmol2_[i]; ii++)
      {
        write_intra(i, ii, density_bins_, natmol2_, intram_mat_density_);
        write_inter_same(i, ii, density_bins_, natmol2_, interm_same_mat_density_, interm_same_maxcdf_mol_);
      }
      for (std::size_t j = i + 1; j < natmol2_.size(); j++)
      {
        for (int ii = 0; ii < natmol2_[i]; ii++)
        {
          write_inter_cross(i, j, ii, density_bins_, natmol2_, cross_index_, interm_cross_mat_density_, interm_cross_maxcdf_mol_);
        }
      }
    }
  }
  printf("finished\n");
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
