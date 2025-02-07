#ifndef _CMDATA_DENSIY_HPP
#define _CMDATA_DENSIY_HPP

#include <cmath>
#include <mutex>
#include <vector>

#include "indexing.hpp"

namespace cmdata::density
{

void kernel_density_estimator(std::vector<float>::iterator x, const std::vector<float> &bins, const double mu, const double norm)
{
  float h = 0.01;
  float from_x = std::max(mu - 2 * h, bins[0]);
  float to_x = std::min(mu + 2 * h, bins.back());
  auto is_geq_start = [&from_x](float i) { return i >= from_x; };
  auto is_geq_end = [&to_x](float i) { return i > to_x; };
  auto start = std::find_if(bins.begin(), bins.end(), is_geq_start);
  auto end = std::find_if(bins.begin(), bins.end(), is_geq_end);
  int from = std::distance(bins.begin(), start);
  int to = std::distance(bins.begin(), end);
  float scale = norm / (0.73853587 * h * std::sqrt(2. * M_PI));
  if (mu < h) scale *= 2.;
  float shift = std::exp(-2.);
  for (int i = from; i < to; i++)
  {
    float f = (mu - bins[i]) / h;
    float kernel = std::exp(-0.5 * f * f);
    x[i] += scale * (kernel - shift);
  }
}

void intra_mol_routine( 
  int i, std::size_t a_i, std::size_t a_j, float dx2, float weight, const std::vector<int> &mol_id_,
  const std::vector<int> &natmol2_, const std::vector<float> &density_bins_,
  const std::vector<float> &inv_num_mol_, std::vector<std::vector<std::mutex>> &frame_same_mutex_, 
  std::vector<std::vector<std::vector<std::vector<float>>>> &intram_mat_density_
)
{
  std::size_t same_mutex_index = cmdata::indexing::mutex_access(mol_id_[i], a_i, a_j, natmol2_);
  std::unique_lock lock(frame_same_mutex_[mol_id_[i]][same_mutex_index]);
  kernel_density_estimator(std::begin(intram_mat_density_[mol_id_[i]][a_i][a_j]), density_bins_, std::sqrt(dx2), weight*inv_num_mol_[i]);
  lock.unlock();
}

void inter_mol_same_routine(
  int i, std::size_t mol_i, std::size_t a_i, std::size_t a_j, float dx2, float weight, 
  const std::vector<int> &mol_id_, const std::vector<int> &natmol2_, const std::vector<float> &density_bins_,
  std::vector<std::vector<std::mutex>> &frame_same_mutex_, std::vector<std::vector<float>> &frame_same_mat_,
  std::vector<std::vector<std::vector<std::vector<float>>>> &interm_same_mat_density_
)
{
  float dist = std::sqrt(dx2);
  std::size_t same_access_index = cmdata::indexing::offset_same(mol_id_[i], mol_i, a_i, a_j, natmol2_);
  std::size_t same_mutex_index = cmdata::indexing::mutex_access(mol_id_[i], a_i, a_j, natmol2_);
  std::unique_lock lock(frame_same_mutex_[mol_id_[i]][same_mutex_index]);
  kernel_density_estimator(std::begin(interm_same_mat_density_[mol_id_[i]][a_i][a_j]), density_bins_, dist, weight);
  frame_same_mat_[mol_id_[i]][same_access_index] = std::min(dist, frame_same_mat_[mol_id_[i]][same_access_index]);
  lock.unlock();
}

void inter_mol_cross_routine(
  int i, int j, std::size_t mol_i, std::size_t mol_j, std::size_t a_i, std::size_t a_j, float dx2, float weight,
  const std::vector<int> &mol_id_, const std::vector<int> &natmol2_, const std::vector<std::vector<int>> &cross_index_,
  const std::vector<float> &density_bins_, std::vector<std::vector<std::mutex>> &frame_cross_mutex_,
  std::vector<std::vector<float>> &frame_cross_mat_, std::vector<std::vector<std::vector<std::vector<float>>>> &interm_cross_mat_density_
)
{
  float dist = std::sqrt(dx2);
  std::size_t cross_access_index = cmdata::indexing::offset_cross(mol_id_[i], mol_id_[j], mol_i, mol_j, a_i, a_j, natmol2_);
  std::size_t cross_mutex_index = cmdata::indexing::mutex_access(mol_id_[j], a_i, a_j, natmol2_);
  std::unique_lock lock(frame_cross_mutex_[cross_index_[mol_id_[i]][mol_id_[j]]][cross_mutex_index]);
  kernel_density_estimator(std::begin(interm_cross_mat_density_[cross_index_[mol_id_[i]][mol_id_[j]]][a_i][a_j]), density_bins_, dist, weight);
  frame_cross_mat_[cross_index_[mol_id_[i]][mol_id_[j]]][cross_access_index] = std::min(dist, frame_cross_mat_[cross_index_[mol_id_[i]][mol_id_[j]]][cross_access_index]);
  lock.unlock();
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

#endif // _CMDATA_DENSIY_HPP
