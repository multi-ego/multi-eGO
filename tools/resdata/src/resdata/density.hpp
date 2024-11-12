#ifndef _RESDATA_DENSIY_HPP
#define _RESDATA_DENSIY_HPP

#include <cmath>
#include <mutex>
#include <vector>

#include "indexing.hpp"

namespace resdata::density
{
//RESIDUE STUFF

void AccumulateIntra( 
  const int i,const std::vector<int> &mol_id_,
  std::vector<std::vector<std::vector<double>>> &resmat_intra_d_, std::vector<std::vector<std::vector<double>>> &resmat_intra_p_,
  std::vector<std::vector<std::vector<std::vector<double>>>> &appo_intra_d, std::vector<std::vector<std::vector<std::vector<double>>>> &appo_intra_p,
  const std::vector<int> &num_mol_unique_
)
{
  for(int im = 0; im<num_mol_unique_[i]; im++)
  {
    for(int res_i=0; res_i<resmat_intra_d_[i].size(); res_i++)
    {
      for(int res_j=0; res_j<resmat_intra_d_[i][0].size(); res_j++)
      {

        if(appo_intra_p[i][im][res_i][res_j]>0.)
        {
          resmat_intra_d_[i][res_i][res_j]+=appo_intra_d[i][im][res_i][res_j];
          resmat_intra_p_[i][res_i][res_j]+=1.;
        }
        appo_intra_d[i][im][res_i][res_j] = 100.;
        appo_intra_p[i][im][res_i][res_j] = 0.;
      }
    }
  }
}

void AccumulateInterSame( 
  const int i,const std::vector<int> &mol_id_,
  std::vector<std::vector<std::vector<double>>> &resmat_same_d_, std::vector<std::vector<std::vector<double>>> &resmat_same_p_,
  std::vector<std::vector<std::vector<std::vector<double>>>> &appo_same_d, std::vector<std::vector<std::vector<std::vector<double>>>> &appo_same_p,
  const std::vector<int> &num_mol_unique_

)
{
  for(int im = 0; im<num_mol_unique_[i]; im++)
  {
    for(int res_i=0; res_i<resmat_same_d_[i].size(); res_i++)
    {
      for(int res_j=0; res_j<resmat_same_d_[i][0].size(); res_j++)
      {
        if(appo_same_p[i][im][res_i][res_j]>0.)
        {
          resmat_same_d_[i][res_i][res_j]+=appo_same_d[i][im][res_i][res_j];
          resmat_same_p_[i][res_i][res_j]+=1.;
        }
        appo_same_d[i][im][res_i][res_j] = 100.;
        appo_same_p[i][im][res_i][res_j] = 0.;
      }
    }
  }
}

void AccumulateInterCross( 
  const int i,const std::vector<int> &mol_id_, const std::vector<std::vector<int>> &cross_index_, const std::vector<int> &natmol2_,
  std::vector<std::vector<std::vector<double>>> &resmat_cross_d_, std::vector<std::vector<std::vector<double>>> &resmat_cross_p_,
  std::vector<std::vector<std::vector<std::vector<double>>>> &appo_cross_d, std::vector<std::vector<std::vector<std::vector<double>>>> &appo_cross_p,
  const std::vector<int> &num_mol_unique_

)
{
  for(int jk=i+1; jk < natmol2_.size() ; jk++)
  {
    for(int res_ik=0; res_ik < resmat_cross_p_[cross_index_[i][jk]].size() ; res_ik++)
    {
      for(int res_jk=0; res_jk < resmat_cross_p_[cross_index_[i][jk]][0].size() ; res_jk++)
      {
        double d_appo=100.;
        double count_appo=0.;
        for(int im = 0; im<num_mol_unique_[i]; im++)
        {
          if(appo_cross_p[cross_index_[i][jk]][im][res_ik][res_jk]>0.)
          {
            d_appo=std::min(appo_cross_d[cross_index_[i][jk]][im][res_ik][res_jk], d_appo);
            count_appo+=1;
          }
          appo_cross_d[cross_index_[i][jk]][im][res_ik][res_jk] = 100.;
          appo_cross_p[cross_index_[i][jk]][im][res_ik][res_jk] = 0.;
        }
        if(count_appo>0){
          resmat_cross_d_[cross_index_[i][jk]][res_ik][res_jk] += d_appo;
          resmat_cross_p_[cross_index_[i][jk]][res_ik][res_jk] += 1.;
        }
      }
    }
  }
}

void NormilizeInterSame( 
  const int n_x_,
  std::vector<std::vector<std::vector<double>>> &resmat_same_d_, std::vector<std::vector<std::vector<double>>> &resmat_same_p_, const std::vector<int> &num_mol_unique_
)
{
  for(int i=0; i<resmat_same_d_.size(); i++)
    {
      for(int res_i=0; res_i<resmat_same_d_[i].size(); res_i++)
      {
        for(int res_j=0; res_j<resmat_same_d_[i][0].size(); res_j++)
        {
          if(resmat_same_p_[i][res_i][res_j]>0.)
          {
            resmat_same_d_[i][res_i][res_j] = resmat_same_d_[i][res_i][res_j]/resmat_same_p_[i][res_i][res_j];///num_mol_unique_[i];
            resmat_same_p_[i][res_i][res_j] = resmat_same_p_[i][res_i][res_j]/n_x_/num_mol_unique_[i];
          }
          else
          {
          resmat_same_d_[i][res_i][res_j] = 0.;
          resmat_same_p_[i][res_i][res_j] = 0.;
          }
        }
      }
    }
}

void NormilizeIntra( 
  const int n_x_,
  std::vector<std::vector<std::vector<double>>> &resmat_intra_d_, std::vector<std::vector<std::vector<double>>> &resmat_intra_p_, const std::vector<int> &num_mol_unique_
)
{
  for(int i=0; i<resmat_intra_d_.size(); i++)
  {
    for(int res_i=0; res_i<resmat_intra_d_[i].size(); res_i++)
    {
      for(int res_j=0; res_j<resmat_intra_d_[i][0].size(); res_j++)
      {
        if(resmat_intra_p_[i][res_i][res_j]>0.)
        {
          resmat_intra_d_[i][res_i][res_j] = resmat_intra_d_[i][res_i][res_j]/resmat_intra_p_[i][res_i][res_j];///num_mol_unique_[i];
          resmat_intra_p_[i][res_i][res_j] = resmat_intra_p_[i][res_i][res_j]/n_x_/num_mol_unique_[i];
        }
        else
        {
        resmat_intra_d_[i][res_i][res_j] = 0.;
        resmat_intra_p_[i][res_i][res_j] = 0.;
        }
      }
    }
  }
}

void NormilizeInterCross( 
  const int n_x_,
  std::vector<std::vector<std::vector<double>>> &resmat_cross_d_, std::vector<std::vector<std::vector<double>>> &resmat_cross_p_
)
{
  for(int i=0; i<resmat_cross_d_.size(); i++)
      {
        for(int res_i=0; res_i<resmat_cross_d_[i].size(); res_i++)
        {
          for(int res_j=0; res_j<resmat_cross_d_[i][0].size(); res_j++)
          {
            if(resmat_cross_p_[i][res_i][res_j]>0.)
            {
            resmat_cross_d_[i][res_i][res_j] = resmat_cross_d_[i][res_i][res_j]/resmat_cross_p_[i][res_i][res_j];
            resmat_cross_p_[i][res_i][res_j] = resmat_cross_p_[i][res_i][res_j]/n_x_;
            }
            else
            {
            resmat_cross_d_[i][res_i][res_j] = 0.;
            resmat_cross_p_[i][res_i][res_j] = 0.;
            }
          }
        }
      }
}


} // namespace resdata::density

#endif // _RESDATA_DENSIY_HPP
