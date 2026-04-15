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
  // Canonicalize to upper triangle: lo <= hi
  if (a_i > a_j) { std::size_t tmp = a_i; a_i = a_j; a_j = tmp; }
  const std::size_t N = static_cast<std::size_t>(natmol2[i]);
  // Row-major upper-triangular index: row lo starts at lo*(2N-lo+1)/2
  const std::size_t per_mol = N * (N + 1) / 2;
  return mol_i * per_mol + a_i * (2 * N - a_i + 1) / 2 + (a_j - a_i);
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
