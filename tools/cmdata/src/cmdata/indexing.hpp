#ifndef _CMDATA_INDEXING_HPP
#define _CMDATA_INDEXING_HPP

#include <vector>

namespace cmdata::indexing
{

int n_bins(const float cut, const float factor = 4.0)
{
  return cut / (0.01 / factor);
}

static std::size_t offset_same( std::size_t i, std::size_t mol_i, std::size_t a_i, std::size_t a_j, const std::vector<int> &natmol2 )
{
  return mol_i * (natmol2[i] * natmol2[i]) + a_i * natmol2[i] + a_j;
}

static std::size_t offset_cross(
  std::size_t i, std::size_t j, std::size_t mol_i, std::size_t mol_j,
  std::size_t a_i, std::size_t a_j, const std::vector<int> &natmol2, std::size_t num_mol_j
)
{
  return mol_i * (num_mol_j * natmol2[i] * natmol2[j]) + mol_j * (natmol2[i] * natmol2[j]) + a_i * natmol2[j] + a_j;
}

} // namespace cmdata::indexing

#endif // _CMDATA_INDEXING_HPP
