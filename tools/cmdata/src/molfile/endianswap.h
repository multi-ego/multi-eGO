/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
University of Illinois Open Source License
Copyright 2003 Theoretical and Computational Biophysics Group, 
All rights reserved.

Developed by:		Theoretical and Computational Biophysics Group
			University of Illinois at Urbana-Champaign
			http://www.ks.uiuc.edu/

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the Software), to deal with 
the Software without restriction, including without limitation the rights to 
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
of the Software, and to permit persons to whom the Software is furnished to 
do so, subject to the following conditions:

Redistributions of source code must retain the above copyright notice, 
this list of conditions and the following disclaimers.

Redistributions in binary form must reproduce the above copyright notice, 
this list of conditions and the following disclaimers in the documentation 
and/or other materials provided with the distribution.

Neither the names of Theoretical and Computational Biophysics Group, 
University of Illinois at Urbana-Champaign, nor the names of its contributors 
may be used to endorse or promote products derived from this Software without 
specific prior written permission.

THE SOFTWARE IS PROVIDED AS IS, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL 
THE CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR 
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
OTHER DEALINGS WITH THE SOFTWARE.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifndef __molfile_endianswap_h
#define __molfile_endianswap_h
/***************************************************************************
 *cr
 *cr            (C) Copyright 1995-2016 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/
/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: endianswap.h,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.8 $       $Date: 2020/10/21 18:03:15 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *   Byte swapping routines used in various plugins
 *   There are two versions of each routine, one that's safe to use in
 *   all cases (but is slow) and one that is only safe to use on memory 
 *   addresses that are aligned to the word size that's being byte-swapped
 *   but are much much much faster.  Use the aligned versions of these
 *   routines whenever possible.  The 'ndata' length count parameters and
 *   internal loops should be safe to use on huge memory arrays on 64-bit
 *   machines.
 *
 ***************************************************************************/

#ifndef ENDIAN_SWAP_H
#define ENDIAN_SWAP_H

#include <stddef.h>


/* Only works with aligned 4-byte quantities, will cause a bus error */
/* on some platforms if used on unaligned data.                      */
static void swap4_aligned(void *v, ptrdiff_t ndata) {
  int *data = (int *) v;
  ptrdiff_t i;
  int *N;
  for (i=0; i<ndata; i++) {
    N = data + i;
    *N=(((*N>>24)&0xff) | ((*N&0xff)<<24) | 
        ((*N>>8)&0xff00) | ((*N&0xff00)<<8));
  }
}


/* Only works with aligned 8-byte quantities, will cause a bus error */
/* on some platforms if used on unaligned data.                      */
static void swap8_aligned(void *v, ptrdiff_t ndata) {
  /* Use int* internally to prevent bugs caused by some compilers */
  /* and hardware that would potentially load data into an FP reg */
  /* and hose everything, such as the old "jmemcpy()" bug in NAMD */
  int *data = (int *) v;  
  ptrdiff_t i;
  int *N; 
  int t0, t1;

  for (i=0; i<ndata; i++) {
    N = data + (i<<1);
    t0 = N[0];
    t0=(((t0>>24)&0xff) | ((t0&0xff)<<24) | 
        ((t0>>8)&0xff00) | ((t0&0xff00)<<8));

    t1 = N[1];
    t1=(((t1>>24)&0xff) | ((t1&0xff)<<24) | 
        ((t1>>8)&0xff00) | ((t1&0xff00)<<8));

    N[0] = t1; 
    N[1] = t0; 
  }
}

#endif
#endif
