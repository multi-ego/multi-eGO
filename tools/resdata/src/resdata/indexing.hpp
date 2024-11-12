#ifndef _RESDATA_INDEXING_HPP
#define _RESDATA_INDEXING_HPP

#include <vector>

namespace resdata::indexing
{

std::size_t mutex_access( std::size_t i, std::size_t r_i, std::size_t r_j, const std::vector<int> &num_res_per_molecule )
{
  return r_i * num_res_per_molecule[i] + r_j;
}

} // namespace resdata::indexing

#endif // _RESDATA_indexing_HPP