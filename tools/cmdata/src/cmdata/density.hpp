#ifndef _CMDATA_DENSITY_HPP
#define _CMDATA_DENSITY_HPP

#include <algorithm>
#include <cmath>
#include <vector>

#include "indexing.hpp"

namespace cmdata::density
{

void kernel_density_estimator(std::vector<float>::iterator x, const std::vector<float> &bins, const float mu, const float norm)
{
  // KDE bandwidth (nm). Controls the smoothing width of each sample.
  static constexpr float h = 0.01f;
  // erf(1/sqrt(2)) — normalisation factor for a Gaussian kernel truncated at ±2h.
  static constexpr float ERF_1_OVER_SQRT2 = 0.73853587f;
  // Kernel tail value at the truncation boundary (exp(-0.5 * 2^2) = exp(-2)).
  static const float KDE_TAIL_SHIFT = std::exp(-2.f);

  // Bins are uniformly spaced: bins[k] = bins[0] + k*dx, so the window
  // [mu-2h, mu+2h] maps to integer indices via O(1) arithmetic instead of
  // two O(n_bins) linear scans with find_if.
  //
  // from: first k where bins[k] >= mu-2h  → ceil((mu-2h - bins[0]) / dx)
  // to:   first k where bins[k] >  mu+2h  → floor((mu+2h - bins[0]) / dx) + 1
  // (using ceil for 'from' is essential: floor gives an extra bin to the left
  //  where kernel < KDE_TAIL_SHIFT, producing a spurious negative accumulation)
  const int n = static_cast<int>(bins.size());
  const float dx = (n > 1) ? (bins[1] - bins[0]) : 1.f;
  const float inv_dx = 1.f / dx;
  int from = static_cast<int>(std::ceil( (mu - 2.f * h - bins[0]) * inv_dx));
  int to   = static_cast<int>(           (mu + 2.f * h - bins[0]) * inv_dx) + 1;
  if (from < 0) from = 0;
  if (to > n)   to   = n;

  float scale = norm / (ERF_1_OVER_SQRT2 * h * std::sqrt(2.f * static_cast<float>(M_PI)));
  if (mu < h) scale *= 2.f;
  for (int i = from; i < to; i++)
  {
    float f = (mu - bins[i]) / h;
    float kernel = std::exp(-0.5f * f * f);
    x[i] += scale * (kernel - KDE_TAIL_SHIFT);
  }
}

void intra_mol_routine(
  int i, std::size_t a_i, std::size_t a_j, float dx2, float weight, const std::vector<int> &mol_id_,
  const std::vector<int> &natmol2_, const std::vector<float> &density_bins_,
  const std::vector<float> &inv_num_mol_,
  std::vector<std::vector<std::vector<std::vector<float>>>> &intram_mat_density_
)
{
  kernel_density_estimator(std::begin(intram_mat_density_[mol_id_[i]][a_i][a_j]), density_bins_, std::sqrt(dx2), weight*inv_num_mol_[i]);
}

void inter_mol_same_routine(
  int i, std::size_t mol_i, std::size_t a_i, std::size_t a_j, float dx2, float weight,
  const std::vector<int> &mol_id_, const std::vector<int> &natmol2_, const std::vector<float> &density_bins_,
  std::vector<std::vector<float>> &frame_same_mat_,
  std::vector<std::vector<std::vector<std::vector<float>>>> &interm_same_mat_density_
)
{
  float dist = std::sqrt(dx2);
  std::size_t same_access_index = cmdata::indexing::offset_same(mol_id_[i], mol_i, a_i, a_j, natmol2_);
  kernel_density_estimator(std::begin(interm_same_mat_density_[mol_id_[i]][a_i][a_j]), density_bins_, dist, weight);
  frame_same_mat_[mol_id_[i]][same_access_index] = std::min(dist, frame_same_mat_[mol_id_[i]][same_access_index]);
}

void inter_mol_cross_routine(
  int i, int j, std::size_t mol_i, std::size_t mol_j, std::size_t a_i, std::size_t a_j, float dx2, float weight,
  const std::vector<int> &mol_id_, const std::vector<int> &natmol2_, const std::vector<std::vector<int>> &cross_index_,
  const std::vector<float> &density_bins_, const std::vector<int> &num_mol_unique,
  std::vector<std::vector<float>> &frame_cross_mat_, std::vector<std::vector<std::vector<std::vector<float>>>> &interm_cross_mat_density_
)
{
  float dist = std::sqrt(dx2);
  std::size_t cross_access_index = cmdata::indexing::offset_cross(mol_id_[i], mol_id_[j], mol_i, mol_j, a_i, a_j, natmol2_, num_mol_unique[mol_id_[j]]);
  kernel_density_estimator(std::begin(interm_cross_mat_density_[cross_index_[mol_id_[i]][mol_id_[j]]][a_i][a_j]), density_bins_, dist, weight);
  frame_cross_mat_[cross_index_[mol_id_[i]][mol_id_[j]]][cross_access_index] = std::min(dist, frame_cross_mat_[cross_index_[mol_id_[i]][mol_id_[j]]][cross_access_index]);
}

void normalize_histo(
    std::size_t i, int ii, int jj, float norm, float inv_num_mol_same,
    std::vector<std::vector<std::vector<std::vector<float>>>> &data
)
{
  std::transform(
    std::begin(data[i][ii][jj]), std::end(data[i][ii][jj]), std::begin(data[i][ii][jj]),
    [norm, inv_num_mol_same]( auto &c ){ return c * norm * inv_num_mol_same; }
  );
}

} // namespace cmdata::density

#endif // _CMDATA_DENSITY_HPP
