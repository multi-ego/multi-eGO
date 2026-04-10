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
#include <xdrfile_trr.h>

#include <algorithm>
#include <cctype>
#include <stdexcept>
#include <string>

namespace cmdata::xtc {

// ---------------------------------------------------------------------------
// Supported trajectory formats
// ---------------------------------------------------------------------------
enum class Format { XTC, TRR, GRO, PDB };

// Infer the trajectory format from the file extension (case-insensitive).
// Throws std::runtime_error for unrecognised extensions.
inline Format detect_format(const std::string &path)
{
    auto dot = path.rfind('.');
    if (dot == std::string::npos)
        throw std::runtime_error(
            "Cannot determine trajectory format: no extension in '" + path + "'");
    std::string ext = path.substr(dot + 1);
    std::transform(ext.begin(), ext.end(), ext.begin(),
                   [](unsigned char c){ return std::tolower(c); });
    if (ext == "xtc")                  return Format::XTC;
    if (ext == "trr")                  return Format::TRR;
    if (ext == "gro")                  return Format::GRO;
    if (ext == "pdb" || ext == "ent")  return Format::PDB;
    throw std::runtime_error(
        "Unsupported trajectory format '." + ext + "'. "
        "Supported formats: .xtc, .trr, .gro, .pdb");
}

// ---------------------------------------------------------------------------
// Frame — lightweight data holder for one trajectory frame.
//
// Memory note: `x` and `offsets` are managed by the owning CMData object and
// must NOT be freed in a Frame destructor (the destructor is intentionally
// omitted so that raw malloc/free in CMData remains safe).
// ---------------------------------------------------------------------------
class Frame {
public:
    int           natom   = 0;
    unsigned long nframe  = 0;        // total frames (0 = unknown; XTC only)
    int64_t      *offsets = nullptr;  // seek table (XTC only; nullptr otherwise)
    int           step    = 0;
    float         time    = 0.f;
    matrix        box     = {};
    rvec         *x       = nullptr;
    float         prec    = 0.f;

    Frame()        { x = static_cast<rvec *>(malloc(0));                   }
    Frame(int n) : natom(n)
                 { x = static_cast<rvec *>(malloc(n * sizeof(rvec)));      }
    // ~Frame() intentionally omitted — x/offsets freed by CMData::~CMData.

    // -----------------------------------------------------------------------
    // XTC: read one frame via xdrfile.
    // Returns exdrOK while frames remain, exdrENDOFFILE at end of file.
    // -----------------------------------------------------------------------
    int read_next_frame(XDRFILE *xd, bool nopbc, PbcType pbc_type, t_pbc *pbc)
    {
        int status = read_xtc(xd, natom, &step, &time, box, x, &prec);
        if (status == exdrOK && !nopbc && pbc) set_pbc(pbc, pbc_type, box);
        return status;
    }

    // -----------------------------------------------------------------------
    // TRR: read one frame via xdrfile.
    // Returns exdrOK while frames remain, exdrENDOFFILE at end of file.
    // Note: this fork's read_trr has an extra has_prop output flag compared to
    // the upstream xdrfile API.
    // -----------------------------------------------------------------------
    int read_next_frame_trr(XDRFILE *xd, bool nopbc, PbcType pbc_type, t_pbc *pbc)
    {
        float   lambda   = 0.f;
        uint8_t has_prop = 0;
        int status = read_trr(xd, natom, &step, &time, &lambda, box, x, nullptr, nullptr, &has_prop);
        if (status == exdrOK && !nopbc && pbc) set_pbc(pbc, pbc_type, box);
        return status;
    }
};

} // namespace cmdata::xtc

#endif // _XTC_FRAME_HPP
