#ifndef _CMDATA_INDEXING_HPP
#define _CMDATA_INDEXING_HPP

#include <vector>

namespace cmdata::indexing
{

int n_bins(const float cut, const float factor = 4.0)
{
  return cut / (0.01 / factor);
}

std::size_t mutex_access( std::size_t i, std::size_t a_i, std::size_t a_j, const std::vector<int> &natmol2 )
{
  return a_i * natmol2[i] + a_j;
}

static std::size_t offset_same( std::size_t i, std::size_t mol_i, std::size_t a_i, std::size_t a_j, const std::vector<int> &natmol2 )
{
  return mol_i * (natmol2[i] * natmol2[i]) + a_i * natmol2[i] + a_j;
}

static std::size_t offset_cross( 
  std::size_t i, std::size_t j, std::size_t mol_i, std::size_t mol_j,
  std::size_t a_i, std::size_t a_j, const std::vector<int> &natmol2
)
{
  return mol_i * (mol_j * natmol2[i] * natmol2[j]) + mol_j * (natmol2[i] * natmol2[j]) + a_i * natmol2[j] + a_j; 
}

struct SameThreadIndices {
  std::size_t start_mti_same;
  std::size_t start_im_same;
  std::size_t start_i_same;
  std::size_t start_j_same;
  std::size_t end_mti_same;
  std::size_t end_im_same;
  std::size_t end_i_same;
  std::size_t end_j_same;
  long int n_loop_operations_same;
};

struct CrossThreadIndices {
  std::size_t start_mti_cross;
  std::size_t start_mtj_cross;
  std::size_t start_im_cross;
  std::size_t start_jm_cross;
  std::size_t start_i_cross;
  std::size_t start_j_cross;
  std::size_t end_mti_cross;
  std::size_t end_mtj_cross;
  std::size_t end_im_cross;
  std::size_t end_jm_cross;
  std::size_t end_i_cross;
  std::size_t end_j_cross;
  int n_loop_operations_cross;
};

} // namespace cmdata::indexing

#endif // _CMDATA_indexing_HPP