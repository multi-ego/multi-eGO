#ifndef _RESDATA_IO_HPP
#define _RESDATA_IO_HPP

#include <gromacs/trajectoryanalysis/topologyinformation.h>
#include <gromacs/fileio/tpxio.h>

#include <filesystem>
#include <string>
#include <fstream>
#include <sstream>
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

namespace resdata::io
{

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

void f_write_inter_cross_res(const std::string &output_prefix,
  std::size_t i, std::size_t j, const std::vector<int> &num_res_per_molecule,
  const std::vector<std::vector<int>> &cross_index,
  const std::vector<int> &residue_indeces_,
  const std::vector<std::vector<std::vector<double>>> &resmat_cross_d_,
  const std::vector<std::vector<std::vector<double>>> &resmat_cross_p_
)
{
  std::filesystem::path ffh = output_prefix + "inter_res_mol_" + std::to_string(i + 1) + "_" + std::to_string(j + 1) + ".dat";
  std::ofstream fp(ffh);

  for ( std::size_t res_i = 0; res_i < num_res_per_molecule[i]; res_i++ )
  {
    for ( std::size_t res_j = 0; res_j < num_res_per_molecule[j]; res_j++ )
    {
    fp << std::fixed << std::setprecision(7) << i + 1;
    fp << " " << std::fixed << std::setprecision(7) << res_i + 1;
    fp << " " << std::fixed << std::setprecision(7) << j + 1;
    fp << " " << std::fixed << std::setprecision(7) << res_j + 1;
    fp << " " << std::fixed << std::setprecision(7) << resmat_cross_d_[cross_index[i][j]][res_i][res_j];
    fp << " " << std::fixed << std::setprecision(7) << resmat_cross_p_[cross_index[i][j]][res_i][res_j];
    fp << "\n";
    }
  }
  fp.close();
}

void f_write_inter_cross_res_block(const std::string &output_prefix,
  std::size_t i, std::size_t j, const std::vector<int> &num_res_per_molecule,
  const std::vector<std::vector<int>> &cross_index,
  const std::vector<int> &residue_indeces_,
  const std::vector<std::vector<std::vector<double>>> &resmat_cross_p_,
  const std::vector<std::vector<std::vector<std::vector<double>>>> &block_resmat_cross_std_,
  const std::vector<int> &block_sizes_
)
{
  std::filesystem::path ffh = output_prefix + "inter_res_mol_" + std::to_string(i + 1) + "_" + std::to_string(j + 1) + ".blocks.dat";
  std::ofstream fp(ffh);
  //header
  fp <<"# mol_i  res_i  mol_j  res_j  p_ij  {blocks_p_ij}";
  fp <<"\n#";
  for (std::size_t i_bk=0; i_bk < block_sizes_.size(); i_bk++)
  {
    fp << " " << std::fixed << std::setprecision(1) << block_sizes_[i_bk];
  }
  fp << "\n";
  for ( std::size_t res_i = 0; res_i < num_res_per_molecule[i]; res_i++ )
  {
    for ( std::size_t res_j = 0; res_j < num_res_per_molecule[j]; res_j++ )
    {
    fp << std::fixed << std::setprecision(7) << i + 1;
    fp << " " << std::fixed << std::setprecision(7) << res_i + 1;
    fp << " " << std::fixed << std::setprecision(7) << j + 1;
    fp << " " << std::fixed << std::setprecision(7) << res_j + 1;
    fp << " " << std::fixed << std::setprecision(7) << resmat_cross_p_[cross_index[i][j]][res_i][res_j];
    // fp << " ; ";
    for (std::size_t i_bk=0; i_bk < block_sizes_.size(); i_bk++)
    {
      fp << " " << std::fixed << std::setprecision(7) << block_resmat_cross_std_[cross_index[i][j]][i_bk][res_i][res_j];
    }
    fp << "\n";
    }
  }
  fp.close();
}

void f_write_inter_same_res(const std::string &output_prefix,
  std::size_t i, const std::vector<int> &num_res_per_molecule,
  //const std::vector<std::vector<int>> &mol_id_,
  const std::vector<int> &residue_indeces_,
  const std::vector<std::vector<std::vector<double>>> &resmat_same_d_,
  const std::vector<std::vector<std::vector<double>>> &resmat_same_p_
)
{
  std::filesystem::path ffh = output_prefix + "inter_res_mol_" + std::to_string(i + 1) + "_" + std::to_string(i + 1) + ".dat";
  std::ofstream fp(ffh);

  for ( std::size_t res_i = 0; res_i < num_res_per_molecule[i]; res_i++ )
  {
    for ( std::size_t res_j = 0; res_j < num_res_per_molecule[i]; res_j++ )
    {
    fp << std::fixed << std::setprecision(7) << i + 1;
    fp << " " << std::fixed << std::setprecision(7) << res_i + 1;
    fp << " " << std::fixed << std::setprecision(7) << i + 1;
    fp << " " << std::fixed << std::setprecision(7) << res_j + 1;
    fp << " " << std::fixed << std::setprecision(7) << resmat_same_d_[i][res_i][res_j];
    fp << " " << std::fixed << std::setprecision(7) << resmat_same_p_[i][res_i][res_j];
    fp << "\n";
    }
  }
  fp.close();
}


void f_write_inter_same_res_block(const std::string &output_prefix,
  std::size_t i, const std::vector<int> &num_res_per_molecule,
  const std::vector<int> &residue_indeces_,
  const std::vector<std::vector<std::vector<double>>> &resmat_same_p_,
  const std::vector<std::vector<std::vector<std::vector<double>>>> &block_resmat_same_std_,
  const std::vector<int> &block_sizes_
)
{
  std::filesystem::path ffh = output_prefix + "inter_res_mol_" + std::to_string(i + 1) + "_" + std::to_string(i + 1) + ".blocks.dat";
  std::ofstream fp(ffh);
  //header
  fp <<"# mol_i  res_i  mol_j  res_j  p_ij  {blocks_p_ij}";
  fp <<"\n#";
  for (std::size_t i_bk=0; i_bk < block_sizes_.size(); i_bk++)
  {
    fp << " " << std::fixed << std::setprecision(1) << block_sizes_[i_bk];
  }
  fp << "\n";

  for ( std::size_t res_i = 0; res_i < num_res_per_molecule[i]; res_i++ )
  {
    for ( std::size_t res_j = 0; res_j < num_res_per_molecule[i]; res_j++ )
    {
    fp << std::fixed << std::setprecision(7) << i + 1;
    fp << " " << std::fixed << std::setprecision(7) << res_i + 1;
    fp << " " << std::fixed << std::setprecision(7) << i + 1;
    fp << " " << std::fixed << std::setprecision(7) << res_j + 1;
    fp << " " << std::fixed << std::setprecision(7) << resmat_same_p_[i][res_i][res_j];
    // fp << " ; ";
    for (std::size_t i_bk=0; i_bk < block_sizes_.size(); i_bk++)
    {
      fp << " " << std::fixed << std::setprecision(7) << block_resmat_same_std_[i][i_bk][res_i][res_j];
    }
    fp << "\n";
    }
  }
  fp.close();
}

void f_write_intra_res(const std::string &output_prefix,
  std::size_t i, const std::vector<int> &num_res_per_molecule,
  //const std::vector<std::vector<int>> &mol_id_,
  const std::vector<int> &residue_indeces_,
  const std::vector<std::vector<std::vector<double>>> &resmat_same_d_,
  const std::vector<std::vector<std::vector<double>>> &resmat_same_p_
)
{
  std::filesystem::path ffh = output_prefix + "intra_res_mol_" + std::to_string(i + 1) + "_" + std::to_string(i + 1) + ".dat";
  std::ofstream fp(ffh);

  for ( std::size_t res_i = 0; res_i < num_res_per_molecule[i]; res_i++ )
  {
    for ( std::size_t res_j = 0; res_j < num_res_per_molecule[i]; res_j++ )
    {
    fp << std::fixed << std::setprecision(7) << i + 1;
    fp << " " << std::fixed << std::setprecision(7) << res_i + 1;
    fp << " " << std::fixed << std::setprecision(7) << i + 1;
    fp << " " << std::fixed << std::setprecision(7) << res_j + 1;
    fp << " " << std::fixed << std::setprecision(7) << resmat_same_d_[i][res_i][res_j];
    fp << " " << std::fixed << std::setprecision(7) << resmat_same_p_[i][res_i][res_j];
    fp << "\n";
    }
  }
  fp.close();
}

void f_write_intra_res_block(const std::string &output_prefix,
  std::size_t i, const std::vector<int> &num_res_per_molecule,
  const std::vector<int> &residue_indeces_,
  const std::vector<std::vector<std::vector<double>>> &resmat_intra_p_,
  const std::vector<std::vector<std::vector<std::vector<double>>>> &block_resmat_intra_std_,
  const std::vector<int> &block_sizes_
)
{
  std::filesystem::path ffh = output_prefix + "intra_res_mol_" + std::to_string(i + 1) + "_" + std::to_string(i + 1) + ".blocks.dat";
  std::ofstream fp(ffh);
  //header
  fp <<"# mol_i  res_i  mol_j  res_j  p_ij  {blocks_p_ij}";
  fp <<"\n#";
  for (std::size_t i_bk=0; i_bk < block_sizes_.size(); i_bk++)
  {
    fp << " " << std::fixed << std::setprecision(1) << block_sizes_[i_bk];
  }
  fp << "\n";

  for ( std::size_t res_i = 0; res_i < num_res_per_molecule[i]; res_i++ )
  {
    for ( std::size_t res_j = 0; res_j < num_res_per_molecule[i]; res_j++ )
    {
    fp << std::fixed << std::setprecision(7) << i + 1;
    fp << " " << std::fixed << std::setprecision(7) << res_i + 1;
    fp << " " << std::fixed << std::setprecision(7) << i + 1;
    fp << " " << std::fixed << std::setprecision(7) << res_j + 1;
    fp << " " << std::fixed << std::setprecision(7) << resmat_intra_p_[i][res_i][res_j];
    // fp << " ; ";
    for (std::size_t i_bk=0; i_bk < block_sizes_.size(); i_bk++)
    {
      fp << " " << std::fixed << std::setprecision(7) << block_resmat_intra_std_[i][i_bk][res_i][res_j];
    }
    fp << "\n";
    }
  }
  fp.close();
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

} // namespace resdata::io

#endif // _CMDATA_IO_HPP


void DefineResidueIndexing(int* index_, const gmx_mtop_t  *mtop_,
  std::vector<int> &num_mol_unique_,std::vector<int> &natmol2_,
  const std::vector<int> &num_mol_unique_temp,const std::vector<int> &natmol2_temp,
  std::vector<int> &num_res_per_molecule,std::vector<int> &residue_indeces_ )
{
  // CREATE VECTOR of RESIDUE indeces to map atom to residue
    int res_ii_appo = 0;
    int molb=0;
    int molb_res=0;
    const char * atomname_appo;
    const char * residuename_appo;

    // for (int i = 0; i < natmol2_.size(); i++)
    //    std::cout << natmol2_[i]*num_mol_unique_[i] << " ";
    int count=0;
    mtopGetAtomAndResidueName(*mtop_, index_[0], &molb, &atomname_appo, &res_ii_appo, &residuename_appo, &molb_res);
    int prevres = res_ii_appo;
    int rescount = 1, i_mol=1;


    // RESIDUE FINDING
    for (int i_mol_type = 0; i_mol_type < num_mol_unique_.size(); i_mol_type++)
    {
      if(i_mol_type>0)
      {
        count += num_mol_unique_[i_mol_type-1]*natmol2_[i_mol_type-1];
      }
      for (int i = 0; i < num_mol_unique_[i_mol_type]*natmol2_[i_mol_type]; i++)
      {
        mtopGetAtomAndResidueName(*mtop_, index_[i+count], &molb, &atomname_appo, &res_ii_appo, &residuename_appo, &molb_res);
        // std::cout<<res_ii_appo<<" "<<residuename_appo<<std::endl;
        if(res_ii_appo != prevres)
        {
          if(i==i_mol*natmol2_[i_mol_type])
          {
            i_mol+=1;
            rescount = 1;
            prevres = res_ii_appo;
          }
          else{
          rescount+=1;
          prevres = res_ii_appo;
          }
        }
        // residue_indeces_.push_back(rescount);      
      }
      //printf("  mol id: %i --> num res:%i\n", i_mol_type, rescount);
      num_res_per_molecule.push_back(rescount);      

      rescount = 0;
      i_mol=1;
    }  

    count=0;
    mtopGetAtomAndResidueName(*mtop_, 0, &molb, &atomname_appo, &res_ii_appo, &residuename_appo, &molb_res);
    prevres = res_ii_appo;
    rescount = 1;
    i_mol=1;

    for (int i_mol_type = 0; i_mol_type < num_mol_unique_temp.size(); i_mol_type++)
    {
      if(i_mol_type>0)
      {
        count += num_mol_unique_temp[i_mol_type-1]*natmol2_temp[i_mol_type-1];
      }
      for (int i = 0; i < num_mol_unique_temp[i_mol_type]*natmol2_temp[i_mol_type]; i++)
      {
        mtopGetAtomAndResidueName(*mtop_, i+count, &molb, &atomname_appo, &res_ii_appo, &residuename_appo, &molb_res);

        if(res_ii_appo != prevres)
        {
          if(i==i_mol*natmol2_temp[i_mol_type])
          {
            i_mol+=1;
            rescount = 1;
            prevres = res_ii_appo;
          }
          else{
          rescount+=1;
          prevres = res_ii_appo;
          }
        }
        residue_indeces_.push_back(rescount);      
      }

      rescount = 0;
      i_mol=1;
    }  
}

