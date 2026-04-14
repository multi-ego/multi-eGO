#ifndef _CMDATA_MINDIST_HPP
#define _CMDATA_MINDIST_HPP

#include <cmath>
#include <vector>

#include "cmdata/indexing.hpp"
#include "cmdata/density.hpp"

namespace cmdata::mindist
{

static void mindist_same(
  const std::vector<float> &density_bins,
  const std::vector<int> &num_mol_unique, const std::vector<int> &natmol2,
  const std::vector<std::vector<float>> &frame_same_mat,
  std::vector<std::vector<std::vector<std::vector<float>>>> &interm_same_maxcdf_mol,
  float weight
)
{
  for (std::size_t mt_i = 0; mt_i < natmol2.size(); mt_i++)
    for (std::size_t im = 0; im < static_cast<std::size_t>(num_mol_unique[mt_i]); im++)
      for (std::size_t i = 0; i < static_cast<std::size_t>(natmol2[mt_i]); i++)
        for (std::size_t j = i; j < static_cast<std::size_t>(natmol2[mt_i]); j++)
        {
          std::size_t offset = cmdata::indexing::offset_same(mt_i, im, i, j, natmol2);
          float mindist = frame_same_mat[mt_i][offset];
          cmdata::density::kernel_density_estimator(std::begin(interm_same_maxcdf_mol[mt_i][i][j]), density_bins, mindist, weight);
          interm_same_maxcdf_mol[mt_i][j][i] = interm_same_maxcdf_mol[mt_i][i][j];
        }
}

static void mindist_cross(
  const std::vector<int> &natmol2, const std::vector<std::vector<int>> &cross_index,
  const std::vector<float> &density_bins, const std::vector<int> &num_mol_unique,
  const std::vector<std::vector<float>> &frame_cross_mat,
  std::vector<std::vector<std::vector<std::vector<float>>>> &interm_cross_maxcdf_mol,
  float weight
)
{
  for (std::size_t mt_i = 0; mt_i < natmol2.size(); mt_i++)
    for (std::size_t mt_j = mt_i + 1; mt_j < natmol2.size(); mt_j++)
      for (std::size_t im = 0; im < static_cast<std::size_t>(num_mol_unique[mt_i]); im++)
        for (std::size_t jm = 0; jm < static_cast<std::size_t>(num_mol_unique[mt_j]); jm++)
          for (std::size_t i = 0; i < static_cast<std::size_t>(natmol2[mt_i]); i++)
            for (std::size_t j = 0; j < static_cast<std::size_t>(natmol2[mt_j]); j++)
            {
              std::size_t offset = cmdata::indexing::offset_cross(mt_i, mt_j, im, jm, i, j, natmol2, num_mol_unique[mt_j]);
              float mindist = frame_cross_mat[cross_index[mt_i][mt_j]][offset];
              cmdata::density::kernel_density_estimator(std::begin(interm_cross_maxcdf_mol[cross_index[mt_i][mt_j]][i][j]), density_bins, mindist, weight);
            }
}

} // namespace cmdata::mindist

#endif // _CMDATA_MINDIST_HPP
