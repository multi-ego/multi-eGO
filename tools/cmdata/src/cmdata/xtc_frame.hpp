#ifndef _XTC_FRAME_HPP
#define _XTC_FRAME_HPP

#ifdef GMXVGE2026
#include "gromacs/utility/vec.h"
#include "gromacs/utility/vectypes.h"
#else
#include <gromacs/math/vec.h>
#include "gromacs/math/vectypes.h"
#endif
#include <gromacs/pbcutil/pbc.h>

#include <xdrfile_xtc.h>

namespace cmdata::xtc {

// ---------------------------------------------------------------------------
// Frame — lightweight data holder for one XTC trajectory frame.
//
// Memory note: `x` and `offsets` are managed by the owning CMData object and
// must NOT be freed in a Frame destructor (the destructor is intentionally
// omitted so that raw malloc/free in CMData remains safe).
// ---------------------------------------------------------------------------
class Frame {
public:
  int           natom   = 0;
  unsigned long nframe  = 0;        // total frames (from read_xtc_header)
  int64_t      *offsets = nullptr;  // seek table (from read_xtc_header)
  int           step    = 0;
  float         time    = 0.f;
  matrix        box     = {};
  rvec         *x       = nullptr;
  float         prec    = 0.f;

  Frame()       { x = (rvec*)malloc(0 * sizeof(rvec)); }
  Frame(int n) : natom(n) { x = (rvec*)malloc(n * sizeof(rvec)); }
  // ~Frame() intentionally omitted — x/offsets freed by CMData::~CMData.

  int read_next_frame(XDRFILE *xd, bool nopbc, PbcType pbc_type, t_pbc *pbc)
  {
    int status = read_xtc(xd, natom, &step, &time, box, x, &prec);
    if (status == exdrOK && !nopbc && pbc) set_pbc(pbc, pbc_type, box);
    return status;
  }
};

} // namespace cmdata::xtc

#endif // _XTC_FRAME_HPP
