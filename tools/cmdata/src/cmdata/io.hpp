#ifndef _CMDATA_IO_HPP
#define _CMDATA_IO_HPP

#include <gromacs/trajectoryanalysis/topologyinformation.h>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <vector>
#include <iomanip>
#include <regex>

#define COUT_FLOAT_PREC6 std::fixed << std::setprecision(6)

static inline void mtopGetMolblockIndex(const gmx_mtop_t& mtop,
                                        int               globalAtomIndex,
                                        int*              moleculeBlock,
                                        int*              moleculeIndex,
                                        int*              atomIndexInMolecule)
{
    // GMX_ASSERT(globalAtomIndex >= 0, "The atom index to look up should not be negative");
    // GMX_ASSERT(globalAtomIndex < mtop.natoms, "The atom index to look up should be within range");
    // GMX_ASSERT(moleculeBlock != nullptr, "molBlock can not be NULL");
    // GMX_ASSERT(!mtop.moleculeBlockIndices.empty(), "The moleculeBlockIndices should not be empty");
    // GMX_ASSERT(*moleculeBlock >= 0,
              //  "The starting molecule block index for the search should not be negative");
    // GMX_ASSERT(*moleculeBlock < gmx::ssize(mtop.moleculeBlockIndices),
              //  "The starting molecule block index for the search should be within range");

    /* Search the molecule block index using bisection */
    int molBlock0 = -1;
    int molBlock1 = mtop.molblock.size();

    int globalAtomStart = 0;
    while (TRUE)
    {
        globalAtomStart = mtop.moleculeBlockIndices[*moleculeBlock].globalAtomStart;
        if (globalAtomIndex < globalAtomStart)
        {
            molBlock1 = *moleculeBlock;
        }
        else if (globalAtomIndex >= mtop.moleculeBlockIndices[*moleculeBlock].globalAtomEnd)
        {
            molBlock0 = *moleculeBlock;
        }
        else
        {
            break;
        }
        *moleculeBlock = ((molBlock0 + molBlock1 + 1) >> 1);
    }

    int molIndex = (globalAtomIndex - globalAtomStart)
                   / mtop.moleculeBlockIndices[*moleculeBlock].numAtomsPerMolecule;
    if (moleculeIndex != nullptr)
    {
        *moleculeIndex = molIndex;
    }
    if (atomIndexInMolecule != nullptr)
    {
        *atomIndexInMolecule = globalAtomIndex - globalAtomStart
                               - molIndex * mtop.moleculeBlockIndices[*moleculeBlock].numAtomsPerMolecule;
    }
}
void mtopGetAtomAndResidueName(const gmx_mtop_t& mtop,
                                             int               globalAtomIndex,
                                             int*              moleculeBlock,
                                             const char**      atomName,
                                             int*              residueNumber,
                                             const char**      residueName,
                                             int*              globalResidueIndex)
{
    int moleculeIndex       = 0;
    int atomIndexInMolecule = 0;
    mtopGetMolblockIndex(mtop, globalAtomIndex, moleculeBlock, &moleculeIndex, &atomIndexInMolecule);

    const gmx_molblock_t&       molb    = mtop.molblock[*moleculeBlock];
    const t_atoms&              atoms   = mtop.moltype[molb.type].atoms;
    const MoleculeBlockIndices& indices = mtop.moleculeBlockIndices[*moleculeBlock];
    if (atomName != nullptr)
    {
        *atomName = *(atoms.atomname[atomIndexInMolecule]);
    }
    if (residueNumber != nullptr)
    {
        if (atoms.nres > mtop.maxResiduesPerMoleculeToTriggerRenumber())
        {
            *residueNumber = atoms.resinfo[atoms.atom[atomIndexInMolecule].resind].nr;
        }
        else
        {
            /* Single residue molecule, keep counting */
            *residueNumber = indices.residueNumberStart + moleculeIndex * atoms.nres
                             + atoms.atom[atomIndexInMolecule].resind;
        }
    }
    if (residueName != nullptr)
    {
        *residueName = *(atoms.resinfo[atoms.atom[atomIndexInMolecule].resind].name);
    }
    if (globalResidueIndex != nullptr)
    {
        *globalResidueIndex = indices.globalResidueStart + moleculeIndex * atoms.nres
                              + atoms.atom[atomIndexInMolecule].resind;
    }
}

namespace cmdata::io
{

void read_symmetry_indices(
  const std::string &path,
  const gmx_mtop_t *top,
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
    throw std::runtime_error(errorMessage.c_str());
  }
  
  int molb = 0;
  std::string residue_entry, atom_entry_i, atom_entry_j, left;
  std::string line;

  while (std::getline(infile, line)) 
  {
    std::istringstream iss(line);
    if (!(iss >> residue_entry >> atom_entry_i >> atom_entry_j)) // each necessary field is there
    {
      if (line=="") continue;
      printf("Skipping line\n%s\n due to syntax non-conformity\n", line.c_str());
      continue;
    }
    if((iss >> left)) printf("Found a problem while reading the symmestry file: %s\n This element is ignored.\n Multiple equivament atoms should be set listing the relevant combinations\n", left.c_str());

    const char *atom_name_i, *atom_name_j, *residue_name_i, *residue_name_j;
    int a_i, a_j, resn_i, resn_j;
    for (std::size_t i = 0; i < natmol2_.size(); i++)
    {
      a_i = 0;
      for (int ii = start_index[i]; ii < start_index[i]+natmol2_[i]; ii++)
      {
        mtopGetAtomAndResidueName(*top, ii, &molb, &atom_name_i, &resn_i, &residue_name_i, nullptr);
	      a_j = 0;
        for (int jj = start_index[i]; jj < start_index[i]+natmol2_[i]; jj++)
        {
          mtopGetAtomAndResidueName(*top, jj, &molb, &atom_name_j, &resn_j, &residue_name_j, nullptr);
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


std::vector<double> read_weights_file( const std::string &path )
{
  std::ifstream infile(path);
  if (!infile.good())
  {
    std::string errorMessage = "Cannot find the indicated weights file";
    throw std::runtime_error(errorMessage.c_str());
  }
  std::vector<double> w;

  std::string line;
  while ( std::getline(infile, line) )
  {
    std::string value;
    std::istringstream iss(line);
    if (line == "")
    {
      printf("Detected empty line. Skipping...\n");
      continue;
    }
    iss >> value;
    w.push_back(std::stod(value));
  }

  if (w.size() == 0)
  {
    std::string errorMessage = "The weights file is empty";
    throw std::runtime_error(errorMessage.c_str());
  }

  for ( std::size_t i = 0; i < w.size(); i++ )
  {
    if (w[i] < 0)
    {
      std::string errorMessage = "The weights file contains negative values";
      throw std::runtime_error(errorMessage.c_str());
    }
  }

  return w;
}


// gmx::RangePartitioning gmx_mtop_molecules(const gmx_mtop_t& mtop)
// {
//     gmx::RangePartitioning mols;

//     for (const gmx_molblock_t& molb : mtop.molblock)
//     {
//         int numAtomsPerMolecule = mtop.moltype[molb.type].atoms.nr;
//         for (int mol = 0; mol < molb.nmol; mol++)
//         {
//             mols.appendBlock(numAtomsPerMolecule);
//         }
//     }

//     return mols;
// }

void f_write_intra(const std::string &output_prefix,
  std::size_t i, int ii, const std::vector<double> &density_bins, const std::vector<int> &natmol2,
  const std::vector<std::vector<std::vector<std::vector<double>>>> &intram_mat_density
)
{
  std::filesystem::path ffh_intra = output_prefix + "intra_mol_" + std::to_string(i + 1) + "_" + std::to_string(i + 1) + "_aa_" + std::to_string(ii + 1) + ".dat";
  std::ofstream fp_intra(ffh_intra);
  for ( std::size_t k = 0; k < density_bins.size(); k++ )
  {
    fp_intra << COUT_FLOAT_PREC6 << density_bins[k];
    for (int jj = 0; jj < natmol2[i]; jj++)
    {
      fp_intra << " " << COUT_FLOAT_PREC6 << intram_mat_density[i][ii][jj][k];
    }
    fp_intra << "\n";
  }

  fp_intra.close();
}

void f_write_inter_same(const std::string &output_prefix,
  std::size_t i, int ii, const std::vector<double> &density_bins, const std::vector<int> &natmol2,
  const std::vector<std::vector<std::vector<std::vector<double>>>> &interm_same_mat_density,
  const std::vector<std::vector<std::vector<std::vector<double>>>> &interm_same_maxcdf_mol
)
{
  std::filesystem::path ffh_inter = output_prefix + "inter_mol_" + std::to_string(i + 1) + "_" + std::to_string(i + 1) + "_aa_" + std::to_string(ii + 1) + ".dat";
  std::ofstream fp_inter(ffh_inter);
  std::filesystem::path ffh_inter_cum = output_prefix + "inter_mol_c_" + std::to_string(i + 1) + "_" + std::to_string(i + 1) + "_aa_" + std::to_string(ii + 1) + ".dat";
  std::ofstream fp_inter_cum(ffh_inter_cum);
  for ( std::size_t k = 0; k < density_bins.size(); k++ )
  {
    fp_inter << COUT_FLOAT_PREC6 << density_bins[k];
    fp_inter_cum << COUT_FLOAT_PREC6 << density_bins[k];
    for (int jj = 0; jj < natmol2[i]; jj++)
    {
      fp_inter << " " << COUT_FLOAT_PREC6 << interm_same_mat_density[i][ii][jj][k];
      fp_inter_cum << " " << COUT_FLOAT_PREC6 << interm_same_maxcdf_mol[i][ii][jj][k];
    }
    fp_inter << "\n";
    fp_inter_cum << "\n";
  }
  fp_inter.close();
  fp_inter_cum.close();
}

void f_write_inter_cross(const std::string &output_prefix,
  std::size_t i, std::size_t j, int ii, const std::vector<double> &density_bins, const std::vector<int> &natmol2,
  const std::vector<std::vector<int>> &cross_index,
  const std::vector<std::vector<std::vector<std::vector<double>>>> &interm_cross_mat_density,
  const std::vector<std::vector<std::vector<std::vector<double>>>> &interm_cross_maxcdf_mol
)
{
  std::filesystem::path ffh = output_prefix + "inter_mol_" + std::to_string(i + 1) + "_" + std::to_string(j + 1) + "_aa_" + std::to_string(ii + 1) + ".dat";
  std::ofstream fp(ffh);
  std::filesystem::path ffh_cum = output_prefix + "inter_mol_c_" + std::to_string(i + 1) + "_" + std::to_string(j + 1) + "_aa_" + std::to_string(ii + 1) + ".dat";
  std::ofstream fp_cum(ffh_cum);
  for ( std::size_t k = 0; k < interm_cross_mat_density[cross_index[i][j]][ii][0].size(); k++ )
  {
    fp << std::fixed << std::setprecision(7) << density_bins[k];
    fp_cum << std::fixed << std::setprecision(7) << density_bins[k];
    for (int jj = 0; jj < natmol2[j]; jj++)
    {
      fp << " " << std::fixed << std::setprecision(7) << interm_cross_mat_density[cross_index[i][j]][ii][jj][k];
      fp_cum << " " << std::fixed << std::setprecision(7) << interm_cross_maxcdf_mol[cross_index[i][j]][ii][jj][k];
    }
    fp << "\n";
    fp_cum << "\n";
  }
  fp.close();
  fp_cum.close();
}

std::vector<uint> read_selection( const std::string &path, const std::string &selection_name )
{
  bool found = false, finished = false;
  std::ifstream infile(path);
  if (!infile.good())
  {
    std::string errorMessage = "Cannot find the indicated selection file";
    throw std::runtime_error(errorMessage.c_str());
  }
  std::vector<uint> sel;

  std::string line;
  while ( std::getline(infile, line) )
  {
    std::string value;
    std::istringstream iss(line);
    if (line == "") continue;

    // find if regex matches the line (regex is no semicolon followed by 0 or more spaces followed by [ and 0 or more spaces followed by selection_name followed by 0 or more spaces followed by ])
    std::regex re_found("([^;]+)\\s*\\[\\s*" + selection_name + "\\s*\\]");
    // find if regex matches the line (same as above but any selection name)
    std::regex re_finished("([^;]+)\\s*\\[\\s*.*\\s*\\]");
    if (std::regex_search(line, re_finished) && found) finished = true;
    if (std::regex_search(line, re_found)) found = true;
    if (found && !finished)
    {
      // check if value is a number
      while (iss >> value)
      {
        try
        {
          std::stoi(value);
        }
        catch (const std::invalid_argument& e)
        {
          std::cerr << "Error: " << e.what() << " in line " << line << std::endl;
          std::cerr << "The value " << value << " found in " << path << " is not a number" << std::endl;
        }
        sel.push_back(std::stoi(value));
      }
    }
  }

  return sel;
}

/**
   * @brief Print a progress bar to the standard output
   * 
   * Taken from https://stackoverflow.com/a/36315819
   * 
   * @param percentage
  */
void print_progress_bar(float percentage)
{
  constexpr std::size_t PROGRESS_BAR_LENGTH = 60;
  constexpr char PROGRESS_BAR[] = "############################################################";
  
  int val = (int) (percentage * 100);
  int lpad = (int) (percentage * PROGRESS_BAR_LENGTH);
  int rpad = PROGRESS_BAR_LENGTH - lpad;
  
  printf("\r%3d%% [%.*s%*s]", val, lpad, PROGRESS_BAR, rpad, "");
  fflush(stdout);
}

} // namespace cmdata::io

#endif // _CMDATA_IO_HPP