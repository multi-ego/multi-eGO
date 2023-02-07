/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2007, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017 by the GROMACS development team.
 * Copyright (c) 2018,2019,2020, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#include "gmxpre.h"

#include <cmath>

#include <algorithm>
#include <numeric>

#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/matio.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/gmxana/gstat.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/mtop_lookup.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

static void clust_size(const char*             ndx,
                       const char*             trx,
                       const char*             xpm,
                       const char*             xpmw,
                       const char*             ncl,
                       const char*             acl,
                       const char*             mcl,
                       const char*             histo,
                       const char*             histotime,
                       const char*             clustime,
                       const char*             trmatrix,
                       const char*             kmatrix,
                       const char*             tempf,
                       const char*             mcn,
                       gmx_bool                bMol,
                       gmx_bool                bPBC,
                       const char*             tpr,
                       double                  cut,
                       double                  mol_cut,
                       int                     bOndx,
                       int                     nskip,
                       int                     skip_last_nmol,
                       int                     nlevels,
                       t_rgb                   rmid,
                       t_rgb                   rhi,
                       int                     ndf,
                       const gmx_output_env_t* oenv)
{
    FILE *       fp, *gp, *hp, *tp, *cndx;
    int*         index = nullptr;
    int          nindex, natoms;
    t_trxstatus* status;
    rvec *       x = nullptr, *v = nullptr, *xcm = nullptr, dx;
    t_pbc        pbc;
    gmx_bool     bSame, bTPRwarn = TRUE;
    /* Topology stuff */
    t_trxframe    fr;
    TpxFileHeader tpxh;
    gmx_mtop_t*   mtop    = nullptr;
    PbcType       pbcType = PbcType::Unset;
    int           ii, jj;
    double          temp, tfac;
    /* Cluster size distribution (matrix) */
    double** cs_dist = nullptr;
    double** tr_matrix = nullptr;
    double** rate_matrix = nullptr;
    double* norm_matrix = nullptr;
    bool* norm_done = nullptr;
    double   tf, dx2, cut2, mcut2, *t_x = nullptr, *t_y, cmid, cmax, cav, ekin;
    int    i, j, k, ai, aj, ci, cj, nframe, nclust, n_x, max_size = 0;
    int *  clust_index, *index_size, *index_old_size, *clust_size, *clust_written, max_clust_size, max_clust_ind, nav, nhisto;
    t_rgb  rlo          = { 1.0, 1.0, 1.0 };
    int    frameCounter = 0;
    double   frameTime;

    clear_trxframe(&fr, TRUE);
    auto timeLabel = output_env_get_time_label(oenv);
    tf             = output_env_get_time_factor(oenv);
    fp             = xvgropen(ncl, "Number of clusters", timeLabel, "N", oenv);
    gp             = xvgropen(acl, "Average cluster size", timeLabel, "#molecules", oenv);
    hp             = xvgropen(mcl, "Max cluster size", timeLabel, "#molecules", oenv);
    tp             = xvgropen(tempf, "Temperature of largest cluster", timeLabel, "T (K)", oenv);

    if (!read_first_frame(oenv, &status, trx, &fr, TRX_NEED_X | TRX_READ_V))
    {
        gmx_file(trx);
    }

    natoms = fr.natoms;
    x      = fr.x;

    if (tpr)
    {
        mtop = new gmx_mtop_t;
        tpxh = readTpxHeader(tpr, true);
        if (tpxh.natoms != natoms)
        {
            gmx_fatal(FARGS, "tpr (%d atoms) and trajectory (%d atoms) do not match!", tpxh.natoms, natoms);
        }
        pbcType = read_tpx(tpr, nullptr, nullptr, &natoms, nullptr, nullptr, mtop);
    }
    if (ndf <= -1)
    {
        tfac = 1;
    }
    else
    {
        tfac = ndf / (3.0 * natoms);
    }

    gmx::RangePartitioning mols;
    if (bMol)
    {
        if (ndx)
        {
            printf("Using molecules rather than atoms. Not reading index file %s\n", ndx);
        }
        GMX_RELEASE_ASSERT(mtop != nullptr, "Trying to access mtop->mols from NULL mtop pointer");
        mols = gmx_mtop_molecules(*mtop);

        /* Make dummy index */
        nindex = mols.numBlocks()-skip_last_nmol;
        snew(index, nindex);
        for (i = 0; (i < nindex); i++)
        {
            index[i] = i;
        }
    }
    else
    {
        char* gname;
        rd_index(ndx, 1, &nindex, &index, &gname);
        sfree(gname);
    }

    snew(clust_index, nindex);
    snew(index_size, nindex);
    snew(index_old_size, nindex);
    snew(clust_size, nindex);
    snew(xcm, nindex);
    /* transition matrix */
    snew(tr_matrix, nindex);
    for(i=0;i<nindex;i++) snew(tr_matrix[i], nindex);
    /* rate matrix */
    snew(rate_matrix, nindex);
    for(i=0;i<nindex;i++) snew(rate_matrix[i], nindex);
    /* norm needed to calculate the rate and transition matrices */
    snew(norm_matrix, nindex);
    /* flag to accumulate correctly the norm matrix */
    snew(norm_done, nindex);
    mcut2 = mol_cut*mol_cut;
    cut2   = cut * cut;
    // total number of trajectory frames
    nframe = 0;
    // number of analysed frames
    n_x    = 0;
    snew(t_y, nindex);
    for (i = 0; (i < nindex); i++)
    {
        t_y[i] = i + 1;
    }
    max_clust_size = 1;
    max_clust_ind  = -1;
    int molb       = 0;
    cndx = xvgropen(clustime, "Index of the oligomer to which each monomer belongs", timeLabel, "Monomer index", oenv);
    double frameTimeStep=1.;
    do
    {
        if(nframe==1&&fr.bTime) frameTimeStep=fr.time;
        if ((nskip == 0) || ((nskip > 0) && ((nframe % nskip) == 0)))
        {
            if (bPBC)
            {
                set_pbc(&pbc, pbcType, fr.box);
            }
            max_clust_size = 1;
            max_clust_ind  = -1;

            /* Put all atoms/molecules in their own cluster, with size 1 */
            for (i = 0; (i < nindex); i++)
            {
                /* Cluster index is indexed with atom index number */
                clust_index[i] = i;
                /* Cluster size is indexed with cluster number */
                clust_size[i] = 1;
                /* Initially each molecule belongs to a cluster of size 1 */
                index_size[i] = 1;
                /* Flag to accumulate the norm matrix */
                norm_done[i] = FALSE;
            }
            /* calculate the center of each molecule */
            for (i = 0; (i < nindex); i++)
            {   
                clear_rvec(xcm[i]);
                ai = index[i];
                double tm = 0.;
                for (ii = mols.block(ai).begin(); ii < mols.block(ai).end(); ii++)
                {
                    for (int m = 0; (m < DIM); m++)
                    {
                        xcm[i][m] += x[ii][m];
                    }
                    tm += 1.0; 
                }
                for (int m = 0; (m < DIM); m++)
                {
                    xcm[i][m] /= tm;
                }
            }

            /* Loop over atoms/molecules */
            for (i = 0; (i < nindex); i++)
            {
                ai = index[i];
                ci = clust_index[i];

                /* Loop over atoms/molecules (only half a matrix) */
                for (j = i + 1; (j < nindex); j++)
                {
                    cj = clust_index[j];

                    if (bPBC)
                    {
                        pbc_dx(&pbc, xcm[i], xcm[j], dx);
                    }
                    else
                    {
                        rvec_sub(xcm[i], xcm[j], dx);
                    }
                    dx2   = iprod(dx, dx);

                    if (dx2 > mcut2) continue;

                    /* If they are not in the same cluster already */
                    if (ci != cj)
                    {
                        aj = index[j];

                        /* Compute distance */
                        if (bMol)
                        {
                            GMX_RELEASE_ASSERT(mols.numBlocks() > 0,
                                               "Cannot access index[] from empty mols");
                            bSame = FALSE;
                            for (ii = mols.block(ai).begin(); !bSame && ii < mols.block(ai).end(); ii++)
                            {
                                for (jj = mols.block(aj).begin(); !bSame && jj < mols.block(aj).end(); jj++)
                                {
                                    if (bPBC)
                                    {
                                        pbc_dx(&pbc, x[ii], x[jj], dx);
                                    }
                                    else
                                    {
                                        rvec_sub(x[ii], x[jj], dx);
                                    }
                                    dx2   = iprod(dx, dx);
                                    bSame = (dx2 < cut2);
                                }
                            }
                        }
                        else
                        {
                            if (bPBC)
                            {
                                pbc_dx(&pbc, x[ai], x[aj], dx);
                            }
                            else
                            {
                                rvec_sub(x[ai], x[aj], dx);
                            }
                            dx2   = iprod(dx, dx);
                            bSame = (dx2 < cut2);
                        }
                        /* If distance less than cut-off */
                        if (bSame)
                        {
                            /* Merge clusters: check for all atoms whether they are in
                             * cluster cj and if so, put them in ci
                             */
                            for (k = 0; (k < nindex); k++)
                            {
                                if (clust_index[k] == cj)
                                {
                                    if (clust_size[cj] <= 0)
                                    {
                                        gmx_fatal(FARGS, "negative cluster size %d for element %d",
                                                  clust_size[cj], cj);
                                    }
                                    clust_size[cj]--;
                                    clust_index[k] = ci;
                                    clust_size[ci]++;
                                }
                            }
                        }
                    }
                }
            }
            for (k = 0; (k < nindex); k++)
            {
                 // this tells how large is the cluster to which each molecule belongs
                 index_size[k] = clust_size[clust_index[k]];
            }
            n_x++;
            srenew(t_x, n_x);
            if (fr.bTime)
            {
                frameTime = fr.time;
            }
            else if (fr.bStep)
            {
                frameTime = fr.step;
            }
            else
            {
                frameTime = ++frameCounter;
            }
            t_x[n_x - 1] = frameTime * tf;
            srenew(cs_dist, n_x);
            snew(cs_dist[n_x - 1], nindex);
            nclust = 0;
            cav    = 0;
            nav    = 0;
            for (i = 0; (i < nindex); i++)
            {
                ci = clust_size[i];
                if (ci > max_clust_size)
                {
                    max_clust_size = ci;
                    max_clust_ind  = i;
                }
                if (ci > 0)
                {
                    nclust++;
                    /* this is the cluster size time-resolved distribution 
                       that is cs[frame][i]=# of oligomers of order (i+1) */
                    cs_dist[n_x - 1][ci - 1] += 1.0;
                    max_size = std::max(max_size, ci);
                    if (ci > 1)
                    {
                        cav += ci;
                        nav++;
                    }
                }
            }
            fprintf(fp, "%14.6e  %10d\n", frameTime, nclust);
            if (nav > 0)
            {
                fprintf(gp, "%14.6e  %10.3f\n", frameTime, cav / nav);
            }
            fprintf(hp, "%14.6e  %10d\n", frameTime, max_clust_size);
            /* update the transition matrix */
            if (n_x>1) 
            {
                double Volume = det(fr.box)*0.0006022;  // NA * nm3->m3
                double Volume2 = Volume*Volume;
                for(i=0;i<nindex;i++)
                {
                   // transition from an oligomer of order index_old_size[i] to on of order index_size[i] 
                   if(cs_dist[n_x-2][index_old_size[i]-1]>0.)
                   {
                     tr_matrix[index_size[i]-1][index_old_size[i]-1]+=1./(cs_dist[n_x-2][index_old_size[i]-1]*((double)index_old_size[i]));
                     if(index_old_size[i]>index_size[i]) {
                       /* k_off */
                       /* this is 1/([oligomer]) that are dissociating */
                       rate_matrix[index_size[i]-1][index_old_size[i]-1]+=Volume/(cs_dist[n_x-2][index_old_size[i]-1]*((double)index_old_size[i]));
                     } else if(index_old_size[i]<index_size[i]){
                       /* k_on */
                       /* this is 1/([oligomer_ligand][oligomer_reactants]) that are associating */
                       double fact=0;
                       for(j=0;j<(index_size[i]-index_old_size[i]);j++) {
                         fact+=cs_dist[n_x-2][j]*((double)(j+1));
                       }
                       fact -= index_old_size[i];
                       rate_matrix[index_size[i]-1][index_old_size[i]-1]+=Volume2/(cs_dist[n_x-2][index_old_size[i]-1]*((double)index_old_size[i])*fact);
                     }
 
                     if(!norm_done[index_old_size[i]-1]) norm_matrix[index_old_size[i]-1]+=1.0;
                     norm_done[index_old_size[i]-1] = TRUE;
                   }
                }
            }
            /* save index_size so that it can be used to generate the transition matrix */
            for(i=0;i<nindex;i++) index_old_size[i] = index_size[i]; 
        }
        /* Analyse velocities, if present */
        if (fr.bV)
        {
            if (!tpr)
            {
                if (bTPRwarn)
                {
                    printf("You need a [REF].tpr[ref] file to analyse temperatures\n");
                    bTPRwarn = FALSE;
                }
            }
            else
            {
                v = fr.v;
                /* Loop over clusters and for each cluster compute 1/2 m v^2 */
                if (max_clust_ind >= 0)
                {
                    ekin = 0;
                    for (i = 0; (i < nindex); i++)
                    {
                        if (clust_index[i] == max_clust_ind)
                        {
                            ai     = index[i];
                            double m = mtopGetAtomMass(mtop, ai, &molb);
                            ekin += 0.5 * m * iprod(v[ai], v[ai]);
                        }
                    }
                    temp = (ekin * 2.0) / (3.0 * tfac * max_clust_size * BOLTZ);
                    fprintf(tp, "%10.3f  %10.3f\n", frameTime, temp);
                }
            }
        }
        fprintf(cndx, "%10.3f ", frameTime);
        for (i = 0; (i < nindex); i++) fprintf(cndx, "%i ", clust_index[i]);
        fprintf(cndx, "\n");

        if((bOndx>1) && (bMol)) 
        {
            int largest=bOndx;
            /* index file per frame per size */
            snew(clust_written, nindex);
            for(int oligsize=2;oligsize<=largest;oligsize++) {
                std::string ndx_name = "cs_" +std::to_string(oligsize) + "_" + std::to_string(nframe) + ".ndx";
                fp = gmx_ffopen(ndx_name.c_str(), "w");
                for (int i = 0; (i < nindex); i++)
                {
                    // this tells how large is the cluster to which each molecule belongs
	            ci = clust_index[i];
	            if(clust_written[ci]==1) continue;
                    if(index_size[i] == oligsize) 
                    {
                        fprintf(fp, "[ clust %i ]\n", ci);
                        for (int j : mols.block(i)) fprintf(fp, "%d\n", j+1);
	                for(int j=i+1; (j < nindex); j++)
	                {
	                    if(clust_index[j]==ci) 
	                    {
                                for (int k : mols.block(j)) fprintf(fp, "%d\n", k+1);
                            }
	                }
	                clust_written[ci]=1;
	            }
	        }
                gmx_ffclose(fp);
            }
            sfree(clust_written);
        }

        nframe++;
    } while (read_next_frame(oenv, status, &fr));
    close_trx(status);
    done_frame(&fr);
    xvgrclose(fp);
    xvgrclose(gp);
    xvgrclose(hp);
    xvgrclose(tp);
    xvgrclose(cndx); 

    snew(clust_written, nindex);
    if (max_clust_ind >= 0)
    {
        fp = gmx_ffopen(mcn, "w");
      /* CARLO: this adds the indices for all the clusters at the end of the trajectory */
      if (bMol)
      {
        for (int i = 0; (i < nindex); i++)
        {
	  ci = clust_index[i];
	  if(clust_written[ci]==1) continue;
          fprintf(fp, "[ clust %i ]\n", ci);
          for (int j : mols.block(i))
          {
             fprintf(fp, "%d\n", j+1);
          }
	  for(int j=i+1; (j < nindex); j++)
	  {
	    if(clust_index[j]==ci) 
	    {
              for (int k : mols.block(j))
              {
                fprintf(fp, "%d\n", k+1);
              }
	    }
	  }
	  clust_written[ci]=1;
	}
      }
        fprintf(fp, "[ max_clust ]\n");
        for (i = 0; (i < nindex); i++)
        {
            if (clust_index[i] == max_clust_ind)
            {
                if (bMol)
                {
                    GMX_RELEASE_ASSERT(mols.numBlocks() > 0,
                                       "Cannot access index[] from empty mols");
                    for (int j : mols.block(i))
                    {
                        fprintf(fp, "%d\n", j + 1);
                    }
                }
                else
                {
                    fprintf(fp, "%d\n", index[i] + 1);
                }
            }
        }
        gmx_ffclose(fp);
    }

    /* Print the double distribution cluster-size/numer, averaged over the trajectory. */
    fp     = xvgropen(histo, "Cluster size distribution", "Cluster size", "()", oenv);
    nhisto = 0;
    fprintf(fp, "%5d  %8.3f\n", 0, 0.0);
    for (j = 0; (j < max_size); j++)
    {
        double nelem = 0;
        for (i = 0; (i < n_x); i++)
        {
            nelem += cs_dist[i][j];
        }
        fprintf(fp, "%5d  %8.3f\n", j + 1, nelem / n_x);
        nhisto += static_cast<int>((j + 1) * nelem / n_x);
    }
    fprintf(fp, "%5d  %8.3f\n", j + 1, 0.0);
    xvgrclose(fp);

    fp = xvgropen(histotime, "Time Resolved distribution of oligomers order", timeLabel, "# of oligomers of order #", oenv);
    for (i = 0; (i < n_x); i++)
    {
        fprintf(fp, "%14.6e ", t_x[i]);
    	for (j = 0; (j < max_size); j++)
    	{
        	fprintf(fp, " %8.3f", cs_dist[i][j]);
        }
        fprintf(fp,"\n");
    }
    xvgrclose(fp); 

    fp = xvgropen(trmatrix, "Transition Matrix", "Oligomers order", "Oligomers order", oenv);
    /* The sum of the rows should be divisible for the oligomer order (that is the row number)\n");
       Rows are transitions toward lower order oligomers\n");
       Columns are transitions toward higher order oligomers\n"); */
    for (i = 0; (i < nindex); i++)
    {
    	for (j = 0; (j < nindex); j++)
    	{
        	fprintf(fp, "%8.6lf ", tr_matrix[i][j]/norm_matrix[j]);
        }
        fprintf(fp,"\n");
    }
    xvgrclose(fp);
 
    fp = xvgropen(kmatrix, "Rates Matrix in ps-1", "Oligomers order", "Oligomers order", oenv);
    /* Rows are transitions toward lower order oligomers\n");
       Columns are transitions toward higher order oligomers\n"); */
    for (i = 0; (i < nindex); i++)
    {
    	for (j = 0; (j < nindex); j++)
    	{
        	fprintf(fp, "%8.6lf ", rate_matrix[i][j]/norm_matrix[j]/frameTimeStep);
        }
        fprintf(fp,"\n");
    }
    xvgrclose(fp);

    fprintf(stderr, "Total number of atoms in clusters =  %d\n", nhisto);

    /* Look for the smallest entry that is not zero
     * This will make that zero is white, and not zero is coloured.
     */
    cmid = 100.0;
    cmax = 0.0;
    for (i = 0; (i < n_x); i++)
    {
        for (j = 0; (j < max_size); j++)
        {
            if ((cs_dist[i][j] > 0) && (cs_dist[i][j] < cmid))
            {
                cmid = cs_dist[i][j];
            }
            cmax = std::max(cs_dist[i][j], cmax);
        }
    }
    fprintf(stderr, "cmid: %g, cmax: %g, max_size: %d\n", cmid, cmax, max_size);
    cmid = 1;
    fp   = gmx_ffopen(xpm, "w");
    //write_xpm3(fp, 0, "Cluster size distribution", "# clusters", timeLabel, "Size", n_x, max_size,
    //           (real)t_x, (real)t_y, cs_dist, 0, cmid, cmax, rlo, rmid, rhi, &nlevels);
    gmx_ffclose(fp);
    cmid = 100.0;
    cmax = 0.0;
    for (i = 0; (i < n_x); i++)
    {
        for (j = 0; (j < max_size); j++)
        {
            cs_dist[i][j] *= (j + 1);
            if ((cs_dist[i][j] > 0) && (cs_dist[i][j] < cmid))
            {
                cmid = cs_dist[i][j];
            }
            cmax = std::max(cs_dist[i][j], cmax);
        }
    }
    fprintf(stderr, "cmid: %g, cmax: %g, max_size: %d\n", cmid, cmax, max_size);
    fp = gmx_ffopen(xpmw, "w");
    //write_xpm3(fp, 0, "Weighted cluster size distribution", "Fraction", timeLabel, "Size", n_x,
    //           max_size, t_x, t_y, cs_dist, 0, cmid, cmax, rlo, rmid, rhi, &nlevels);
    gmx_ffclose(fp);
    delete mtop;
    sfree(t_x);
    sfree(t_y);
    for (i = 0; (i < n_x); i++)
    {
        sfree(cs_dist[i]);
    }
    sfree(cs_dist);
    sfree(clust_index);
    sfree(clust_size);
    sfree(index);
}

static void do_interm_mat(const char*             trx,
                          const char*             outfile_inter,
                          const char*             outfile_intra,
                          gmx_bool                bPBC,
                          const char*             tpr,
                          double                  cut,
                          double                  mol_cut,
                          double                  d_pow,
                          int                     nskip,
                          int                     skip_last_nmol,
                          gmx_bool                write_histo,
                          const gmx_output_env_t* oenv)
{
    t_trxframe    fr;
    clear_trxframe(&fr, TRUE);

    t_trxstatus* status;
    if (!read_first_frame(oenv, &status, trx, &fr, TRX_NEED_X | TRX_READ_V))
    {
        gmx_file(trx);
    }

    int natoms = fr.natoms;
    rvec *x = fr.x;

    TpxFileHeader tpxh;
    gmx_mtop_t    *mtop = nullptr;
    PbcType pbcType = PbcType::Unset;
    if (tpr)
    {
        mtop = new gmx_mtop_t;
        tpxh = readTpxHeader(tpr, true);
        if (tpxh.natoms != natoms)
        {
            gmx_fatal(FARGS, "tpr (%d atoms) and trajectory (%d atoms) do not match!", tpxh.natoms, natoms);
        }
        pbcType = read_tpx(tpr, nullptr, nullptr, &natoms, nullptr, nullptr, mtop);
    }

    gmx::RangePartitioning mols;
    GMX_RELEASE_ASSERT(tpr, "Cannot access topology without having read it from TPR");
    mols = gmx_mtop_molecules(*mtop);

    // number of molecules
    int nindex = mols.numBlocks()-skip_last_nmol;
    std::vector<int> num_mol;
    num_mol.push_back(1);
    int num_unique_molecules=0;
    // number of atoms per molecule, assuming them identical when consecutive molecules have the same number of atoms
    std::vector<int> natmol2;
    natmol2.push_back(mols.block(0).end());
    for(unsigned i=1; i<nindex; i++) {
       natmol2.push_back(mols.block(i).end()-mols.block(i-1).end());
       if(natmol2[i]==natmol2[i-1]) num_mol[num_unique_molecules]++;
       else {
         num_mol.push_back(1);
         num_unique_molecules++;
       }
    }
    std::vector<int>::iterator it = std::unique(natmol2.begin(), natmol2.end());  
    natmol2.resize(std::distance(natmol2.begin(),it));

    std::vector<int> start_index;
    std::vector<int> mol_id;
    mol_id.push_back(0);
    std::vector<double> inv_num_mol;
    start_index.push_back(0); 
    num_unique_molecules=0;
    inv_num_mol.push_back(1./(static_cast<double>(num_mol[num_unique_molecules])));

    for(unsigned i=1; i<nindex; i++) {
       if(mols.block(i).end()-mols.block(i-1).end()==natmol2[num_unique_molecules]) {
          start_index.push_back(start_index[i-1]);
       } else {
          start_index.push_back(natmol2[num_unique_molecules]);
          num_unique_molecules++;
       }
       mol_id.push_back(num_unique_molecules);
       inv_num_mol.push_back(1./static_cast<double>(num_mol[num_unique_molecules]));
    }

    printf("number of different molecules %u\n", natmol2.size());
    for(unsigned i=0; i<natmol2.size();i++) printf("mol %u num %u size %u\n", i, num_mol[i], natmol2[i]);
    //for(unsigned i=0; i<nindex;i++) printf("start_idnex %u %u\n", i, start_index[i]);
    //for(unsigned i=0; i<nindex;i++) printf("invnummol %u %lf\n", i, inv_num_mol[i]);

    // matrix atm x atm for probabilities
    std::vector<std::vector<std::vector<double> > > interm_same_mat(natmol2.size());    
    std::vector<std::vector<std::vector<double> > > interm_cross_mat((natmol2.size()*(natmol2.size()-1))/2);    
    std::vector<std::vector<std::vector<double> > > intram_mat(natmol2.size());
    // matrix atm x atm for normalising distances  
    std::vector<std::vector<std::vector<double> > > interm_same_mat_dist_count(natmol2.size());    
    std::vector<std::vector<std::vector<double> > > interm_cross_mat_dist_count((natmol2.size()*(natmol2.size()-1))/2);    
    std::vector<std::vector<std::vector<double> > > intram_mat_dist_count(natmol2.size());    
    // matrix atm x atm for average distances 
    std::vector<std::vector<std::vector<double> > > interm_same_mat_dist(natmol2.size());    
    std::vector<std::vector<std::vector<double> > > interm_cross_mat_dist((natmol2.size()*(natmol2.size()-1))/2);    
    std::vector<std::vector<std::vector<double> > > intram_mat_dist(natmol2.size());    
    // matrix atm x atm for average r12 distances 
    std::vector<std::vector<std::vector<double> > > interm_same_mat_dist12(natmol2.size());    
    std::vector<std::vector<std::vector<double> > > interm_cross_mat_dist12((natmol2.size()*(natmol2.size()-1))/2);    
    std::vector<std::vector<std::vector<double> > > intram_mat_dist12(natmol2.size());    
    // Tensor atm x atm x 110 to accumulate histograms
    std::vector<std::vector<std::vector<std::vector<int> > > > interm_same_mat_histo(natmol2.size());  
    std::vector<std::vector<std::vector<std::vector<int> > > > interm_cross_mat_histo((natmol2.size()*(natmol2.size()-1))/2);    
    std::vector<std::vector<std::vector<std::vector<int> > > > intram_mat_histo(natmol2.size());

    int cross_count=0;
    std::vector<std::vector<int> > cross_index(natmol2.size(), std::vector<int>(natmol2.size(),0));
    for(int i=0; i<natmol2.size();i++) {
      interm_same_mat[i].resize(natmol2[i], std::vector<double>(natmol2[i],0.));
      intram_mat[i].resize(natmol2[i], std::vector<double>(natmol2[i],0.));
      interm_same_mat_dist_count[i].resize(natmol2[i], std::vector<double>(natmol2[i],0.));
      intram_mat_dist_count[i].resize(natmol2[i], std::vector<double>(natmol2[i],0.));
      interm_same_mat_dist[i].resize(natmol2[i], std::vector<double>(natmol2[i],0.));
      intram_mat_dist[i].resize(natmol2[i], std::vector<double>(natmol2[i],0.));
      interm_same_mat_dist12[i].resize(natmol2[i], std::vector<double>(natmol2[i],0.));
      intram_mat_dist12[i].resize(natmol2[i], std::vector<double>(natmol2[i],0.));
      interm_same_mat_histo[i].resize(natmol2[i], std::vector<std::vector<int>>(natmol2[i], std::vector<int>(55,0)));
      intram_mat_histo[i].resize(natmol2[i], std::vector<std::vector<int>>(natmol2[i], std::vector<int>(55,0)));
      for(int j=i+1; j<natmol2.size();j++) {
        interm_cross_mat[cross_count].resize(natmol2[i], std::vector<double>(natmol2[j],0.));
        interm_cross_mat_dist_count[cross_count].resize(natmol2[i], std::vector<double>(natmol2[j],0.));
        interm_cross_mat_dist[cross_count].resize(natmol2[i], std::vector<double>(natmol2[j],0.));
        interm_cross_mat_dist12[cross_count].resize(natmol2[i], std::vector<double>(natmol2[j],0.));
        interm_cross_mat_histo[cross_count].resize(natmol2[i], std::vector<std::vector<int>>(natmol2[j], std::vector<int>(55,0)));
        cross_index[i][j]=cross_count;
        cross_count++;
      }
    } 

    // vector of center of masses
    rvec *xcm = nullptr;
    snew(xcm, nindex);

    double mcut2 = mol_cut*mol_cut;
    double cut2   = cut * cut;
    // total number of trajectory frames
    int nframe = 0;
    // number of analysed frames
    int n_x = 0;

    printf("Ready!\n"); fflush(stdout);

    do
    {
        if ((nskip == 0) || ((nskip > 0) && ((nframe % nskip) == 0)))
        {
            t_pbc pbc;
            if (bPBC) set_pbc(&pbc, pbcType, fr.box);

            /* calculate the center of each molecule */
            for (int i = 0; (i < nindex); i++)
            {   
                clear_rvec(xcm[i]);
                double tm = 0.;
                for (int ii = mols.block(i).begin(); ii < mols.block(i).end(); ii++)
                {
                    for (int m = 0; (m < DIM); m++)
                    {
                        xcm[i][m] += x[ii][m];
                    }
                    tm += 1.0; 
                }
                for (int m = 0; (m < DIM); m++)
                {
                    xcm[i][m] /= tm;
                }
            }
   
            /* Loop over molecules */
            for (int i = 0; i < nindex; i++)
            {
                // Temporary structures for intermediate values
                // this is to set that at least on interaction has been found
                std::vector<std::vector<int> > added(natmol2[mol_id[i]], std::vector<int>(natmol2[mol_id[i]], 0));
                std::vector<std::vector<std::vector<int> > > added_cross((natmol2.size()*(natmol2.size()-1))/2);
                // matrices atm x atm for accumulating distances 
                std::vector<std::vector<double> > interm_same_mat_mdist(natmol2[mol_id[i]], std::vector<double>(natmol2[mol_id[i]], 100.));    
                std::vector<std::vector<double> > interm_same_mat_Mdist12(natmol2[mol_id[i]], std::vector<double>(natmol2[mol_id[i]], 0.));    
                std::vector<std::vector<double> > intram_mat_mdist(natmol2[mol_id[i]], std::vector<double>(natmol2[mol_id[i]], 100.));    
                std::vector<std::vector<double> > intram_mat_Mdist12(natmol2[mol_id[i]], std::vector<double>(natmol2[mol_id[i]], 0.));    
                std::vector<std::vector<std::vector<double> > > interm_cross_mat_mdist((natmol2.size()*(natmol2.size()-1))/2);
                std::vector<std::vector<std::vector<double> > > interm_cross_mat_Mdist12((natmol2.size()*(natmol2.size()-1))/2); 
                for (int j = mol_id[i]+1; j < natmol2.size(); j++) {
                    interm_cross_mat_mdist[cross_index[mol_id[i]][j]].resize(natmol2[mol_id[i]], std::vector<double>(natmol2[mol_id[j]], 100.));    
                    interm_cross_mat_Mdist12[cross_index[mol_id[i]][j]].resize(natmol2[mol_id[i]], std::vector<double>(natmol2[mol_id[j]], 0.));
                    added_cross[cross_index[mol_id[i]][j]].resize(natmol2[mol_id[i]], std::vector<int>(natmol2[mol_id[j]], 0));
                }
                /* Loop over molecules  */
                for (int j = 0; j < nindex; j++)
                {
                    rvec dx;
                    if (bPBC) pbc_dx(&pbc, xcm[i], xcm[j], dx);
                    else rvec_sub(xcm[i], xcm[j], dx);
                    double dx2 = iprod(dx, dx);
                    if (dx2 > mcut2) continue;
                    if(mol_id[i]!=mol_id[j]&&j<i) continue;

                    /* Compute distance */
                    int a_i = 0;
                    double norm = std::max(inv_num_mol[i],inv_num_mol[j]);
                    GMX_RELEASE_ASSERT(mols.numBlocks() > 0,"Cannot access index[] from empty mols");
                    for (int ii = mols.block(i).begin(); ii < mols.block(i).end(); ii++)
                    {
                        int a_j = 0;
                        for (int jj = mols.block(j).begin(); jj < mols.block(j).end(); jj++)
                        {
                            if (bPBC) pbc_dx(&pbc, x[ii], x[jj], dx);
                            else rvec_sub(x[ii], x[jj], dx);
                            double dx2 = iprod(dx, dx);
                            if(dx2 < cut2) {
                                double id12 = std::pow(1./dx2,-0.5*d_pow);
                                if(i!=j) { // intermolecular
                                   if(mol_id[i]==mol_id[j]) { // inter same molecule specie
                                      if(!added[a_i][a_j]) {
                                        interm_same_mat[mol_id[i]][a_i][a_j] +=  norm;
                                        if(a_i!=a_j) interm_same_mat[mol_id[i]][a_j][a_i] += norm;
                                        added[a_i][a_j] = 1;
                                        added[a_j][a_i] = 1;
                                      }
                                      interm_same_mat_mdist[a_i][a_j] = std::min(interm_same_mat_mdist[a_i][a_j], sqrt(dx2));
                                      interm_same_mat_mdist[a_j][a_i] = std::min(interm_same_mat_mdist[a_i][a_j], sqrt(dx2));
                                      interm_same_mat_Mdist12[a_i][a_j] = std::max(interm_same_mat_Mdist12[a_i][a_j], id12);
                                      interm_same_mat_Mdist12[a_j][a_i] = std::max(interm_same_mat_Mdist12[a_i][a_j], id12);
                                   } else { // inter cross molecule specie
                                      if(!added_cross[cross_index[mol_id[i]][mol_id[j]]][a_i][a_j]) {
                                        interm_cross_mat[cross_index[mol_id[i]][mol_id[j]]][a_i][a_j] +=  norm;
                                        added_cross[cross_index[mol_id[i]][mol_id[j]]][a_i][a_j] = 1;
                                      }
                                      interm_cross_mat_mdist[cross_index[mol_id[i]][mol_id[j]]][a_i][a_j] = std::min(interm_cross_mat_mdist[cross_index[mol_id[i]][mol_id[j]]][a_i][a_j], sqrt(dx2));
                                      interm_cross_mat_Mdist12[cross_index[mol_id[i]][mol_id[j]]][a_i][a_j] = std::max(interm_cross_mat_Mdist12[cross_index[mol_id[i]][mol_id[j]]][a_i][a_j], id12);
                                   }
                                } else { // intramolecular
                                   intram_mat[mol_id[i]][a_i][a_j] += norm;
                                   intram_mat_mdist[a_i][a_j] = std::min(intram_mat_mdist[a_i][a_j], sqrt(dx2));
                                   intram_mat_mdist[a_j][a_i] = std::min(intram_mat_mdist[a_i][a_j], sqrt(dx2));
                                   intram_mat_Mdist12[a_i][a_j] = std::max(intram_mat_Mdist12[a_i][a_j], id12);
                                   intram_mat_Mdist12[a_j][a_i] = std::max(intram_mat_Mdist12[a_i][a_j], id12);
                                }
                            }
                            a_j++;
                        }
                        a_i++;
                    }
                }
                for(int ii=0; ii<natmol2[mol_id[i]]; ii++) {
                   for(int jj=0; jj<natmol2[mol_id[i]]; jj++) {
                      if(interm_same_mat_mdist[ii][jj]<100.) {
                        if(write_histo) interm_same_mat_histo[mol_id[i]][ii][jj][static_cast<unsigned>(std::floor(interm_same_mat_mdist[ii][jj]/(cut/55.)))]++;
                        interm_same_mat_dist[mol_id[i]][ii][jj] += interm_same_mat_mdist[ii][jj];
                        interm_same_mat_dist12[mol_id[i]][ii][jj] += interm_same_mat_Mdist12[ii][jj];
                        interm_same_mat_dist_count[mol_id[i]][ii][jj]+=1.;
                      } 
                      if(intram_mat_mdist[ii][jj]<100.) {
                        if(write_histo) intram_mat_histo[mol_id[i]][ii][jj][static_cast<unsigned>(std::floor(intram_mat_mdist[ii][jj]/(cut/55.)))]++;
                        intram_mat_dist[mol_id[i]][ii][jj] += intram_mat_mdist[ii][jj];
                        intram_mat_dist12[mol_id[i]][ii][jj] += intram_mat_Mdist12[ii][jj];
                        intram_mat_dist_count[mol_id[i]][ii][jj]+=1.;
                      }
                   }
                }
                for (int j = mol_id[i]+1; j < natmol2.size(); j++) {
                   for(int ii=0; ii<natmol2[mol_id[i]]; ii++) {
                      for(int jj=0; jj<natmol2[mol_id[j]]; jj++) {
                         if(interm_cross_mat_mdist[cross_index[mol_id[i]][j]][ii][jj]<100.) {
                           if(write_histo) interm_cross_mat_histo[cross_index[mol_id[i]][j]][ii][jj][static_cast<unsigned>(std::floor(interm_cross_mat_mdist[cross_index[mol_id[i]][j]][ii][jj]/(cut/55.)))]++;
                           interm_cross_mat_dist[cross_index[mol_id[i]][j]][ii][jj] += interm_cross_mat_mdist[cross_index[mol_id[i]][j]][ii][jj];
                           interm_cross_mat_dist12[cross_index[mol_id[i]][j]][ii][jj] += interm_cross_mat_Mdist12[cross_index[mol_id[i]][j]][ii][jj];
                           interm_cross_mat_dist_count[cross_index[mol_id[i]][j]][ii][jj]+=1.;
                         }
                      }
                   }
                } 
            }
            n_x++;
        }
        nframe++;
    } while (read_next_frame(oenv, status, &fr));
    close_trx(status);
    done_frame(&fr);

    sfree(xcm);
    printf("Done!\n"); fflush(stdout);

    // normalisations
    for(int i=0; i<natmol2.size(); i++) {
       for(int ii=0; ii<natmol2[i]; ii++) {
          for(int jj=0; jj<natmol2[i]; jj++) {
             interm_same_mat[i][ii][jj] /= static_cast<double>(n_x);
             intram_mat[i][ii][jj] /= static_cast<double>(n_x);
             if(interm_same_mat_dist_count[i][ii][jj] > 0) {
                interm_same_mat_dist[i][ii][jj] = interm_same_mat_dist[i][ii][jj]/interm_same_mat_dist_count[i][ii][jj];
                interm_same_mat_dist12[i][ii][jj] = std::pow(interm_same_mat_dist12[i][ii][jj]/interm_same_mat_dist_count[i][ii][jj], 1./d_pow);
             } else {
                interm_same_mat_dist[i][ii][jj] = 0.;
                interm_same_mat_dist12[i][ii][jj] = 0.;
             }
             if(intram_mat_dist_count[i][ii][jj] > 0) {
                intram_mat_dist[i][ii][jj] = intram_mat_dist[i][ii][jj]/intram_mat_dist_count[i][ii][jj];
                intram_mat_dist12[i][ii][jj] = std::pow(intram_mat_dist12[i][ii][jj]/intram_mat_dist_count[i][ii][jj], 1./d_pow);
             } else {
                intram_mat_dist[i][ii][jj] = 0.;
                intram_mat_dist12[i][ii][jj] = 0.;
             }
             if(write_histo) {
                FILE *fp = nullptr;
                std::string ffh = "inter_mol_"+std::to_string(i+1)+"_"+std::to_string(i+1)+"_aa_"+std::to_string(ii+1)+"_"+std::to_string(jj+1)+".dat";
                fp = gmx_ffopen(ffh, "w");
                for(int k=0; k<55; k++) {
                   fprintf(fp, "%lf %d\n",  cut/55.*k+cut/110., interm_same_mat_histo[i][ii][jj][k]);
                }
                gmx_ffclose(fp);
                ffh = "intra_mol_"+std::to_string(i+1)+"_"+std::to_string(i+1)+"_aa_"+std::to_string(ii+1)+"_"+std::to_string(jj+1)+".dat";
                fp = gmx_ffopen(ffh, "w");
                for(int k=0; k<55; k++) {
                   fprintf(fp, "%lf %d\n",  cut/55.*k+cut/110., intram_mat_histo[i][ii][jj][k]);
                }
                gmx_ffclose(fp);
             }
          }
       }
       for(int j=i+1; j<natmol2.size(); j++) {
          for(int ii=0; ii<natmol2[i]; ii++) {
             for(int jj=0; jj<natmol2[j]; jj++) {
                interm_cross_mat[cross_index[i][j]][ii][jj] /= static_cast<double>(n_x);
                if(interm_cross_mat_dist_count[cross_index[i][j]][ii][jj] > 0) {
                   interm_cross_mat_dist[cross_index[i][j]][ii][jj] = interm_cross_mat_dist[cross_index[i][j]][ii][jj]/interm_cross_mat_dist_count[cross_index[i][j]][ii][jj];
                   interm_cross_mat_dist12[cross_index[i][j]][ii][jj] = std::pow(interm_cross_mat_dist12[cross_index[i][j]][ii][jj]/interm_cross_mat_dist_count[cross_index[i][j]][ii][jj], 1./d_pow);
                } else {
                   interm_cross_mat_dist[cross_index[i][j]][ii][jj] = 0.;
                   interm_cross_mat_dist12[cross_index[i][j]][ii][jj] = 0.;
                }
                if(write_histo) {
                   FILE *fp = nullptr;
                   std::string ffh = "inter_mol_"+std::to_string(i+1)+"_"+std::to_string(j+1)+"_aa_"+std::to_string(ii+1)+"_"+std::to_string(jj+1)+".dat";
                   fp = gmx_ffopen(ffh, "w");
                   for(int k=0; k<55; k++) {
                      fprintf(fp, "%lf %d\n",  cut/55.*k+cut/110., interm_cross_mat_histo[cross_index[i][j]][ii][jj][k]);
                   }
                   gmx_ffclose(fp);
                }
             }
          }
       }
    }

    for(int i=0; i<natmol2.size(); i++) {
       FILE *fp = nullptr;
       std::string inter_file_name(outfile_inter);
       std::size_t found = inter_file_name.find_last_of(".");
       fp = gmx_ffopen(inter_file_name.insert(found,"_"+std::to_string(i+1)+"_"+std::to_string(i+1)), "w");
       for(int ii=0; ii<natmol2[i]; ii++) {
          for(int jj=0; jj<natmol2[i]; jj++) {
             fprintf(fp, "%4i %4i %4i %4i %9.6lf %9.6lf %9.6lf\n", i+1, ii+1, i+1, jj+1, interm_same_mat_dist[i][ii][jj], interm_same_mat_dist12[i][ii][jj], interm_same_mat[i][ii][jj]);
          }
       }
       gmx_ffclose(fp);
       std::string intra_file_name(outfile_intra);
       found = intra_file_name.find_last_of(".");
       fp = gmx_ffopen(intra_file_name.insert(found,"_"+std::to_string(i+1)+"_"+std::to_string(i+1)), "w");
       for(int ii=0; ii<natmol2[i]; ii++) {
          for(int jj=0; jj<natmol2[i]; jj++) {
             fprintf(fp, "%4i %4i %4i %4i %9.6lf %9.6lf %9.6lf\n", i+1, ii+1, i+1, jj+1, intram_mat_dist[i][ii][jj], intram_mat_dist12[i][ii][jj], intram_mat[i][ii][jj]);
          }
       }
       gmx_ffclose(fp);
       for(int j=i+1; j<natmol2.size(); j++) {
          std::string inter_c_file_name(outfile_inter);
          found = inter_c_file_name.find_last_of(".");
          fp = gmx_ffopen(inter_c_file_name.insert(found,"_"+std::to_string(i+1)+"_"+std::to_string(j+1)), "w");
          for(int ii=0; ii<natmol2[i]; ii++) {
             for(int jj=0; jj<natmol2[j]; jj++) {
                fprintf(fp, "%4i %4i %4i %4i %9.6lf %9.6lf %9.6lf\n", i+1, ii+1, j+1, jj+1, interm_cross_mat_dist[cross_index[i][j]][ii][jj], interm_cross_mat_dist12[cross_index[i][j]][ii][jj], interm_cross_mat[cross_index[i][j]][ii][jj]);
             }
          }
          gmx_ffclose(fp);
       }
    }
}

int gmx_clustsize(int argc, char* argv[])
{
    const char* desc[] = {
        "[THISMODULE] computes the size distributions of molecular/atomic clusters in",
        "the gas phase. The output is given in the form of an [REF].xpm[ref] file.",
        "The total number of clusters is written to an [REF].xvg[ref] file.[PAR]",
        "When the [TT]-mol[tt] option is given clusters will be made out of",
        "molecules rather than atoms, which allows clustering of large molecules.",
        "In this case an index file would still contain atom numbers",
        "or your calculation will die with a SEGV.[PAR]",
        "When velocities are present in your trajectory, the temperature of",
        "the largest cluster will be printed in a separate [REF].xvg[ref] file assuming",
        "that the particles are free to move. If you are using constraints,",
        "please correct the temperature. For instance water simulated with SHAKE",
        "or SETTLE will yield a temperature that is 1.5 times too low. You can",
        "compensate for this with the [TT]-ndf[tt] option. Remember to take the removal",
        "of center of mass motion into account.[PAR]",
        "The [TT]-mc[tt] option will produce an index file containing the",
        "atom numbers of the largest cluster."
    };

    real     cutoff = 0.50;
    real     mol_cutoff = 6.00;
    real     d_pow = -6.;
    int      bOndx   = 0;
    int      nskip   = 0;
    int      skip_last_nmol = 0;
    int      nlevels = 20;
    int      ndf     = -1;
    gmx_bool bMol    = FALSE;
    gmx_bool bPBC    = TRUE;
    gmx_bool iMAT    = FALSE;
    gmx_bool iMAThis = FALSE;
    rvec     rlo     = { 1.0, 1.0, 0.0 };
    rvec     rhi     = { 0.0, 0.0, 1.0 };

    gmx_output_env_t* oenv;

    t_pargs pa[] = {
        { "-cut",
          FALSE,
          etREAL,
          { &cutoff },
          "Largest distance (nm) to be considered in a cluster" },
        { "-mol_cut",
          FALSE,
          etREAL,
          { &mol_cutoff },
          "Largest distance (nm) to be considered between molecules in a cluster" },
        { "-mol",
          FALSE,
          etBOOL,
          { &bMol },
          "Cluster molecules rather than atoms (needs [REF].tpr[ref] file)" },
        { "-inter_mol",
          FALSE,
          etBOOL,
          { &iMAT },
          "Perform an inter/intra-molecular interactions analysis (needs [REF].tpr[ref] file)" },
        { "-histo",
          FALSE,
          etBOOL,
          { &iMAThis },
          "with -inter_mol plots the histogram of the distances for all pairs (needs [REF].tpr[ref] file)" },
        { "-skip_last_nmol", FALSE, etINT, { &skip_last_nmol }, "Number of molecules to skip (from the end)" },
        { "-tr_olig_ndx",
          FALSE,
          etINT,
          { &bOndx },
          "write index files for all oligomers size from 2 to tr_olig_ndx for every frame, it could enerate A LOT of files" },
        { "-d_pow", FALSE, etREAL, { &d_pow }, "Averaging of distance^power" },
        { "-pbc", FALSE, etBOOL, { &bPBC }, "Use periodic boundary conditions" },
        { "-nskip", FALSE, etINT, { &nskip }, "Number of frames to skip between writing" },
        { "-nlevels",
          FALSE,
          etINT,
          { &nlevels },
          "Number of levels of grey in [REF].xpm[ref] output" },
        { "-ndf",
          FALSE,
          etINT,
          { &ndf },
          "Number of degrees of freedom of the entire system for temperature calculation. "
          "If not set, the number of atoms times three is used." },
        { "-rgblo",
          FALSE,
          etRVEC,
          { rlo },
          "RGB values for the color of the lowest occupied cluster size" },
        { "-rgbhi",
          FALSE,
          etRVEC,
          { rhi },
          "RGB values for the color of the highest occupied cluster size" }
    };
#define NPA asize(pa)
    const char *fnNDX, *fnTPR;
    t_rgb       rgblo, rgbhi;

    t_filenm fnm[] = {
        { efTRX, "-f", nullptr, ffREAD },         { efTPR, nullptr, nullptr, ffOPTRD },
        { efNDX, nullptr, nullptr, ffOPTRD },     { efXPM, "-o", "csize", ffWRITE },
        { efXPM, "-ow", "csizew", ffWRITE },      { efXVG, "-nc", "nclust", ffWRITE },
        { efXVG, "-mc", "maxclust", ffWRITE },    { efXVG, "-ac", "avclust", ffWRITE },
        { efXVG, "-hc", "histo-clust", ffWRITE }, { efXVG, "-temp", "temp", ffOPTWR },
        { efXVG, "-hct", "histo-time", ffWRITE },
        { efXVG, "-ict", "clust-index-time", ffWRITE },
        { efXVG, "-trm", "transitions-matrix", ffWRITE },
        { efXVG, "-km", "rates-matrix", ffWRITE },
        { efNDX, "-mcn", "maxclust", ffOPTWR },
        { efNDX, "-irmat", "intermat", ffOPTWR },
        { efNDX, "-iamat", "intramat", ffOPTWR }
    };
#define NFILE asize(fnm)

    if (!parse_common_args(&argc, argv, PCA_CAN_VIEW | PCA_CAN_TIME | PCA_TIME_UNIT, NFILE, fnm,
                           NPA, pa, asize(desc), desc, 0, nullptr, &oenv))
    {
        return 0;
    }

    if(iMAT) bMol = TRUE;
 
    fnNDX   = ftp2fn_null(efNDX, NFILE, fnm);
    rgblo.r = rlo[XX];
    rgblo.g = rlo[YY];
    rgblo.b = rlo[ZZ];
    rgbhi.r = rhi[XX];
    rgbhi.g = rhi[YY];
    rgbhi.b = rhi[ZZ];

    fnTPR = ftp2fn_null(efTPR, NFILE, fnm);
    if (bMol && !fnTPR)
    {
        gmx_fatal(FARGS, "You need a tpr file for the -mol option");
    }

    if(!iMAT)
    clust_size(fnNDX, ftp2fn(efTRX, NFILE, fnm), opt2fn("-o", NFILE, fnm), opt2fn("-ow", NFILE, fnm),
               opt2fn("-nc", NFILE, fnm), opt2fn("-ac", NFILE, fnm), opt2fn("-mc", NFILE, fnm),
               opt2fn("-hc", NFILE, fnm), opt2fn("-hct", NFILE, fnm), opt2fn("-ict", NFILE, fnm), opt2fn("-trm", NFILE, fnm), opt2fn("-km", NFILE, fnm), 
               opt2fn("-temp", NFILE, fnm), opt2fn("-mcn", NFILE, fnm), bMol, bPBC, fnTPR, cutoff, mol_cutoff, bOndx, nskip, skip_last_nmol, nlevels, rgblo, rgbhi, ndf, oenv);

    else if(iMAT)
    do_interm_mat(ftp2fn(efTRX, NFILE, fnm),
                  opt2fn("-irmat", NFILE, fnm),
                  opt2fn("-iamat", NFILE, fnm),
                  bPBC,
                  fnTPR,
                  cutoff,
                  mol_cutoff,
                  d_pow,
                  nskip,
                  skip_last_nmol,
                  iMAThis,
                  oenv);

    output_env_done(oenv);

    return 0;
}
