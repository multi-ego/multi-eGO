#ifndef _RESDATA_BLOCKAVG_HPP
#define _RESDATA_BLOCKAVG_HPP

#include <cmath>
#include <mutex>
#include <vector>

#include "indexing.hpp"

namespace resdata::blockavg
{

void AccumulateInterCrossBlock( 
  const int i,const int n_x_,const std::vector<int> &mol_id_, const std::vector<std::vector<int>> &cross_index_, const std::vector<int> &natmol2_,
  std::vector<std::vector<std::vector<std::vector<double>>>> &block_resmat_cross_p_,std::vector<std::vector<std::vector<std::vector<double>>>> &block_resmat_cross_p2_,std::vector<std::vector<std::vector<std::vector<double>>>> &block_resmat_cross_p_temp, std::vector<std::vector<std::vector<std::vector<int>>>> &block_counter_cross,
  const std::vector<int> &block_sizes_,std::vector<std::vector<std::vector<std::vector<double>>>> &appo_cross_p,
  const std::vector<int> &num_mol_unique_
)
{
  for(int jk=i+1; jk < natmol2_.size() ; jk++)
  {
    for(int i_bk = 0; i_bk<block_resmat_cross_p_[cross_index_[i][jk]].size(); i_bk ++)
    {
      for(int res_ik=0; res_ik < block_resmat_cross_p_[cross_index_[i][jk]][0].size() ; res_ik++)
      {
        for(int res_jk=0; res_jk < block_resmat_cross_p_[cross_index_[i][jk]][0][0].size() ; res_jk++)
        {
          double appo=0.;
          for(int im=0; im<num_mol_unique_[i]; im++)
          {
            if(appo_cross_p[cross_index_[i][jk]][im][res_ik][res_jk]>0.)
            {
              appo+=1;
            }
          }
          if(n_x_>0 && appo>0)block_resmat_cross_p_temp[cross_index_[i][jk]][i_bk][res_ik][res_jk] += 1.;
          if(std::fmod(n_x_, block_sizes_[i_bk]) == 0 && n_x_ > 0)
          {
            double a = block_resmat_cross_p_temp[cross_index_[i][jk]][i_bk][res_ik][res_jk]/block_sizes_[i_bk];
            
            block_resmat_cross_p_[cross_index_[i][jk]][i_bk][res_ik][res_jk] += a;
            block_resmat_cross_p2_[cross_index_[i][jk]][i_bk][res_ik][res_jk] += a*a;

            block_resmat_cross_p_temp[cross_index_[i][jk]][i_bk][res_ik][res_jk] = 0.;
            block_counter_cross[cross_index_[i][jk]][i_bk][res_ik][res_jk]+=1;
          }
        }
      }  
    }
  }
}


void AccumulateInterSameBlock( 
  const int i,const int n_x_,const std::vector<int> &mol_id_, const std::vector<int> &natmol2_,
  std::vector<std::vector<std::vector<std::vector<double>>>> &block_resmat_same_p_,std::vector<std::vector<std::vector<std::vector<double>>>> &block_resmat_same_p2_,std::vector<std::vector<std::vector<std::vector<double>>>> &block_resmat_same_p_temp, std::vector<std::vector<std::vector<std::vector<int>>>> &block_counter_same,
  const std::vector<int> &block_sizes_,std::vector<std::vector<std::vector<std::vector<double>>>> &appo_same_p,
  const std::vector<int> &num_mol_unique_
)
{
  for(int i_bk = 0; i_bk<block_resmat_same_p_[i].size(); i_bk ++)
  {
    for(int res_i=0; res_i<block_resmat_same_p_[i][i_bk].size(); res_i++)
    {
      for(int res_j=0; res_j<block_resmat_same_p_[i][i_bk][0].size(); res_j++)
      {
        for(int im=0; im<num_mol_unique_[i]; im++)
        {
          if(appo_same_p[i][im][res_i][res_j]>0.)
          {
            if(n_x_>0)block_resmat_same_p_temp[i][i_bk][res_i][res_j] += appo_same_p[i][im][res_i][res_j];
            if(std::fmod(n_x_, block_sizes_[i_bk]) == 0 && n_x_ > 0 && i+1==num_mol_unique_[i])
            {
              double a = block_resmat_same_p_temp[i][i_bk][res_i][res_j] / (block_sizes_[i_bk]*num_mol_unique_[i]);
              
              block_resmat_same_p_[i][i_bk][res_i][res_j] += a;
              block_resmat_same_p2_[i][i_bk][res_i][res_j] += a*a;

              block_resmat_same_p_temp[i][i_bk][res_i][res_j] = 0.;
              block_counter_same[i][i_bk][res_i][res_j]+=1;
            }
          }
        }  
      }
    }
  }
}

void AccumulateIntraBlock( 
  const int i,const int n_x_,const std::vector<int> &mol_id_, const std::vector<int> &natmol2_,
  std::vector<std::vector<std::vector<std::vector<double>>>> &block_resmat_intra_p_,std::vector<std::vector<std::vector<std::vector<double>>>> &block_resmat_intra_p2_,std::vector<std::vector<std::vector<std::vector<double>>>> &block_resmat_intra_p_temp, std::vector<std::vector<std::vector<std::vector<int>>>> &block_counter_intra,
  const std::vector<int> &block_sizes_,std::vector<std::vector<std::vector<std::vector<double>>>> &appo_intra_p,
  const std::vector<int> &num_mol_unique_

)
{
  for(int i_bk = 0; i_bk<block_resmat_intra_p_[i].size(); i_bk ++)
  {
    for(int res_i=0; res_i<block_resmat_intra_p_[i][i_bk].size(); res_i++)
    {
      for(int res_j=0; res_j<block_resmat_intra_p_[i][i_bk][0].size(); res_j++)
      {
        for(int im=0; im<num_mol_unique_[i]; im++)
        {
          if(appo_intra_p[i][im][res_i][res_j]>0.)
          {
            if(n_x_>0)block_resmat_intra_p_temp[i][i_bk][res_i][res_j] += appo_intra_p[i][im][res_i][res_j];
            if(std::fmod(n_x_, block_sizes_[i_bk]) == 0 && n_x_ > 0 && i+1==num_mol_unique_[i])
            {
              double a = block_resmat_intra_p_temp[i][i_bk][res_i][res_j] / (block_sizes_[i_bk]*num_mol_unique_[i]);
              
              block_resmat_intra_p_[i][i_bk][res_i][res_j] += a;
              block_resmat_intra_p2_[i][i_bk][res_i][res_j] += a*a;

              block_resmat_intra_p_temp[i][i_bk][res_i][res_j] = 0.;
              block_counter_intra[i][i_bk][res_i][res_j]+=1;
            }
          }
        }
      }
    }
  }
}




void NormilizeInterCrossBlock( 
  std::vector<std::vector<std::vector<std::vector<double>>>> &block_resmat_cross_std_, std::vector<std::vector<std::vector<std::vector<double>>>> &block_resmat_cross_p_,
  std::vector<std::vector<std::vector<std::vector<double>>>> &block_resmat_cross_p2_,std::vector<std::vector<std::vector<std::vector<int>>>> &block_counter_cross
)
{
  for(int i=0; i<block_resmat_cross_p_.size(); i++)
  {
    for(int i_bk=0; i_bk<block_resmat_cross_p_[i].size(); i_bk++)
    {
      for(int res_i=0; res_i<block_resmat_cross_p_[i][i_bk].size(); res_i++)
      {
        for(int res_j=0; res_j<block_resmat_cross_p_[i][i_bk][0].size(); res_j++)
        {
          block_resmat_cross_p_[i][i_bk][res_i][res_j]/=block_counter_cross[i][i_bk][res_i][res_j];
          block_resmat_cross_p2_[i][i_bk][res_i][res_j]/=block_counter_cross[i][i_bk][res_i][res_j]; 
          block_resmat_cross_std_[i][i_bk][res_i][res_j]=std::sqrt((block_resmat_cross_p2_[i][i_bk][res_i][res_j]-block_resmat_cross_p_[i][i_bk][res_i][res_j]*block_resmat_cross_p_[i][i_bk][res_i][res_j])/block_counter_cross[i][i_bk][res_i][res_j]);
        }
      }
    }
  }
}

void NormilizeInterSameBlock( 
  std::vector<std::vector<std::vector<std::vector<double>>>> &block_resmat_same_std_, std::vector<std::vector<std::vector<std::vector<double>>>> &block_resmat_same_p_,
  std::vector<std::vector<std::vector<std::vector<double>>>> &block_resmat_same_p2_,std::vector<std::vector<std::vector<std::vector<int>>>> &block_counter_same,
  const std::vector<int> &num_mol_unique_
)
{
  for(int i=0; i<block_resmat_same_p_.size(); i++)
  {
    for(int i_bk=0; i_bk<block_resmat_same_p_[i].size(); i_bk++)
    {
      for(int res_i=0; res_i<block_resmat_same_p_[i][i_bk].size(); res_i++)
      {
        for(int res_j=0; res_j<block_resmat_same_p_[i][i_bk][0].size(); res_j++)
        {
          if(block_counter_same[i][i_bk][res_i][res_j]>0)
          {
            block_resmat_same_p_[i][i_bk][res_i][res_j]/=block_counter_same[i][i_bk][res_i][res_j];
            block_resmat_same_p2_[i][i_bk][res_i][res_j]/=block_counter_same[i][i_bk][res_i][res_j];
            block_resmat_same_std_[i][i_bk][res_i][res_j]=std::sqrt((block_resmat_same_p2_[i][i_bk][res_i][res_j]-block_resmat_same_p_[i][i_bk][res_i][res_j]*block_resmat_same_p_[i][i_bk][res_i][res_j])/block_counter_same[i][i_bk][res_i][res_j]);
          }
          else block_resmat_same_std_[i][i_bk][res_i][res_j]=0.;

        }
      }
    }
  }
}

void NormilizeIntraBlock( 
  std::vector<std::vector<std::vector<std::vector<double>>>> &block_resmat_intra_std_, std::vector<std::vector<std::vector<std::vector<double>>>> &block_resmat_intra_p_,
  std::vector<std::vector<std::vector<std::vector<double>>>> &block_resmat_intra_p2_,std::vector<std::vector<std::vector<std::vector<int>>>> &block_counter_intra,
  const std::vector<int> &num_mol_unique_
)
{
  for(int i=0; i<block_resmat_intra_p_.size(); i++)
  {
    for(int i_bk=0; i_bk<block_resmat_intra_p_[i].size(); i_bk++)
    {
      for(int res_i=0; res_i<block_resmat_intra_p_[i][i_bk].size(); res_i++)
      {
        for(int res_j=0; res_j<block_resmat_intra_p_[i][i_bk][0].size(); res_j++)
        {
          if(block_counter_intra[i][i_bk][res_i][res_j]>0)
          {
            block_resmat_intra_p_[i][i_bk][res_i][res_j]/=block_counter_intra[i][i_bk][res_i][res_j];
            block_resmat_intra_p2_[i][i_bk][res_i][res_j]/=block_counter_intra[i][i_bk][res_i][res_j];
            block_resmat_intra_std_[i][i_bk][res_i][res_j]=std::sqrt((block_resmat_intra_p2_[i][i_bk][res_i][res_j]-block_resmat_intra_p_[i][i_bk][res_i][res_j]*block_resmat_intra_p_[i][i_bk][res_i][res_j])/block_counter_intra[i][i_bk][res_i][res_j]);
          }
          else block_resmat_intra_std_[i][i_bk][res_i][res_j]=0.;
        
        }
        
      }
    }
  }
}


} // namespace resdata::block_avg

#endif // _RESDATA_BLOCKAVG_HPP
