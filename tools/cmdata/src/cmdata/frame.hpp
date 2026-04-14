#ifndef _CMDATA_FRAME_HPP
#define _CMDATA_FRAME_HPP

#include <cmath>
#include <cstring>
#include <string>
#include <algorithm>
#include <stdexcept>

#ifdef GMXVGE2026
#include "gromacs/utility/vec.h"
#include "gromacs/utility/vectypes.h"
#else
#include <gromacs/math/vec.h>
#include "gromacs/math/vectypes.h"
#endif
#include <gromacs/pbcutil/pbc.h>

#include "molfile/molfile_plugin.h"

namespace cmdata {
namespace traj {

// Status codes returned by Frame::read_next_frame
static constexpr int FRAME_OK  = 0;
static constexpr int FRAME_EOF = 1;

// ---------------------------------------------------------------------------
// Frame — one trajectory frame, backed by a VMD molfile plugin.
//
// coord_scale controls the unit conversion applied to molfile coordinates:
//   1.0f  — plugin already outputs nm (GROMACS formats compiled with MOLFILE_NATIVE_NM)
//   0.1f  — plugin outputs Angstrom (PDB)
//
// The same scale is applied to box lengths (ts.A/B/C); angles are dimensionless.
// ---------------------------------------------------------------------------
class Frame {
public:
  int    natom      = 0;
  float  time       = 0.f;
  float  coord_scale = 0.1f; // Angstrom→nm by default; 1.0f when plugin outputs nm
  matrix box        = {};
  rvec  *x          = nullptr;

  molfile_plugin_t *mf_plugin = nullptr;
  void             *mf_handle = nullptr;
  float            *mf_coords = nullptr; // flat xyzxyz..., in plugin units

  Frame() { x = (rvec *)malloc(0); }

  Frame(int n, molfile_plugin_t *plugin, void *handle, float scale = 0.1f)
    : natom(n), coord_scale(scale), mf_plugin(plugin), mf_handle(handle)
  {
    x         = (rvec *)malloc(n * sizeof(rvec));
    mf_coords = (float *)malloc(3 * n * sizeof(float));
  }

  // Convert crystallographic box (lengths in plugin units; angles in degrees)
  // to GROMACS triclinic matrix (nm).  scale converts lengths to nm.
  static void box_from_molfile(float A, float B, float C,
                               float alpha, float beta, float gamma,
                               float scale, matrix out)
  {
    constexpr float d2r = static_cast<float>(M_PI) / 180.f;
    const float a  = A * scale, b = B * scale, c = C * scale;
    const float ca = std::cos(alpha * d2r);
    const float cb = std::cos(beta  * d2r);
    const float cg = std::cos(gamma * d2r);
    const float sg = std::sin(gamma * d2r);

    out[0][0] = a;        out[0][1] = 0.f;                     out[0][2] = 0.f;
    out[1][0] = b * cg;   out[1][1] = b * sg;                  out[1][2] = 0.f;
    out[2][0] = c * cb;   out[2][1] = c * (ca - cb * cg) / sg;
    const float z2 = c * c - out[2][0] * out[2][0] - out[2][1] * out[2][1];
    out[2][2] = z2 > 0.f ? std::sqrt(z2) : 0.f;
  }

  int read_next_frame(bool nopbc, PbcType pbc_type, t_pbc *pbc)
  {
    molfile_timestep_t ts{};
    ts.coords = mf_coords;
    if (mf_plugin->read_next_timestep(mf_handle, natom, &ts) != MOLFILE_SUCCESS)
      return FRAME_EOF;

    for (int i = 0; i < natom; i++)
    {
      x[i][XX] = mf_coords[3 * i    ] * coord_scale;
      x[i][YY] = mf_coords[3 * i + 1] * coord_scale;
      x[i][ZZ] = mf_coords[3 * i + 2] * coord_scale;
    }
    time = static_cast<float>(ts.physical_time);
    box_from_molfile(ts.A, ts.B, ts.C, ts.alpha, ts.beta, ts.gamma, coord_scale, box);
    if (!nopbc && pbc) set_pbc(pbc, pbc_type, box);
    return FRAME_OK;
  }
};

} // namespace traj
} // namespace cmdata

#endif // _CMDATA_FRAME_HPP
