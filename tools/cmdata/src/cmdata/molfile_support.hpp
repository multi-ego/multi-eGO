#ifndef _CMDATA_MOLFILE_SUPPORT_HPP
#define _CMDATA_MOLFILE_SUPPORT_HPP

#include <cstring>
#include <stdexcept>
#include <string>

#include "molfile/vmdplugin.h"
#include "molfile/molfile_plugin.h"

// ---------------------------------------------------------------------------
// Declarations for the plugin init/register functions after macro renaming.
//
// Each .cpp is compiled with per-source COMPILE_DEFINITIONS:
//   gromacsplugin.cpp: -DVMDPLUGIN_init=molfile_gromacsplugin_init
//                      -DVMDPLUGIN_register=molfile_gromacsplugin_register
//   pdbplugin.cpp:     -DVMDPLUGIN_init=molfile_pdbplugin_init
//                      -DVMDPLUGIN_register=molfile_pdbplugin_register
// ---------------------------------------------------------------------------
int molfile_gromacsplugin_init();
int molfile_gromacsplugin_register(void *, vmdplugin_register_cb);
int molfile_pdbplugin_init();
int molfile_pdbplugin_register(void *, vmdplugin_register_cb);

namespace cmdata {

// Return the molfile plugin whose filename_extension matches `ext`
// (e.g. "gro", "trr", "pdb").  Initialises all plugins on first call.
inline molfile_plugin_t *get_molfile_plugin(const std::string &ext)
{
  static bool inited = false;
  if (!inited)
  {
    molfile_gromacsplugin_init();
    molfile_pdbplugin_init();
    inited = true;
  }

  struct Ctx { const char *ext; molfile_plugin_t *plugin; };
  Ctx ctx{ ext.c_str(), nullptr };

  // Non-capturing lambda → compatible with vmdplugin_register_cb.
  // filename_extension may be a comma-separated list (e.g. "pdb,ent"),
  // so tokenise and compare each token individually.
  auto cb = [](void *data, vmdplugin_t *p) -> int {
    auto *c  = static_cast<Ctx *>(data);
    auto *mp = reinterpret_cast<molfile_plugin_t *>(p);
    if (!mp->filename_extension) return 0;
    // duplicate so we can strtok in-place without mutating plugin data
    char *exts = strdup(mp->filename_extension);
    char *tok  = strtok(exts, ",");
    while (tok)
    {
      if (std::strcmp(c->ext, tok) == 0) { c->plugin = mp; break; }
      tok = strtok(nullptr, ",");
    }
    free(exts);
    return 0;
  };

  molfile_gromacsplugin_register(&ctx, cb);
  molfile_pdbplugin_register(&ctx, cb);

  if (!ctx.plugin)
    throw std::runtime_error("No molfile plugin found for extension: " + ext);
  return ctx.plugin;
}

} // namespace cmdata

#endif // _CMDATA_MOLFILE_SUPPORT_HPP
