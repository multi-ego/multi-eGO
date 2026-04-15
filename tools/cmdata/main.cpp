// standard library imports
#include <filesystem>
#include <iostream>
#include <string>
// CLI11 argument parser (header-only)
#include <CLI11.hpp>
// cmdata import
#include "src/cmdata/cmdata.hpp"

int main(int argc, const char** argv)
{
  std::cout << "################################################" << std::endl;
  std::cout << "#################### cmdata ####################" << std::endl;
  std::cout << "################### Multi-eGO ##################" << std::endl;
  std::cout << "################## Version 1.0 #################" << std::endl;
  std::cout << "################################################\n" << std::endl;

  CLI::App app{"cmdata — contact data from GROMACS trajectories for multi-eGO"};
  app.set_version_flag("--version", "1.0");

  std::string traj_path, top_path;
  std::string mode         = "intra+same+cross";
  std::string out_prefix   = "";
  std::string bkbn_H       = "";
  std::string weights_path = "";
  float  cutoff     = 0.75f;
  float  mol_cutoff = 6.0f;
  float  t_begin    = 0.0f;
  float  t_end      = -1.0f;
  int    nskip      = 0;
  int    dt         = 0;
  bool   nopbc      = false;
  bool   noh5       = false;
  bool   h5         = false;
#ifdef USE_HDF5
  h5 = true;
#endif

  app.add_option("-f,--traj",       traj_path,     "Input trajectory file (.xtc, .trr, .gro, .pdb)")->required()->check(CLI::ExistingFile);
  app.add_option("-s,--top",        top_path,      "Input topology file (.tpr, .pdb, or .gro)")->required()->check(CLI::ExistingFile);
  app.add_option("--mode",          mode,          "Calculation mode: +-separated combination of intra, same, cross");
  app.add_option("-o,--out",        out_prefix,    "Output file prefix");
  app.add_option("-b,--t_begin",    t_begin,       "Start time (ps)");
  app.add_option("-e,--t_end",      t_end,         "End time (ps); -1 means read to the end");
  app.add_option("--dt",            dt,            "Only process frames at multiples of this time (ps)");
  app.add_option("--cutoff",        cutoff,        "Distance cutoff for atom pairs (nm)");
  app.add_option("--mol_cutoff",    mol_cutoff,    "Centre-of-mass cutoff for molecule pairs (nm)");
  app.add_option("--nskip",         nskip,         "Skip every N frames (0 = no skipping)");
  app.add_option("--bkbn_H",        bkbn_H,        "Backbone hydrogen atom name to include (H and HN are always included)");
  app.add_option("--weights",       weights_path,  "Per-frame weight file (one float per line)")->check(CLI::ExistingFile);
  app.add_flag("--no_pbc",          nopbc,         "Disable periodic boundary corrections");
  app.add_flag("--noh5",            noh5,          "Write plain-text .dat files instead of HDF5 .h5");

  CLI11_PARSE(app, argc, argv);

  // --noh5 overrides the HDF5 default set at compile time
  if (noh5) h5 = false;

  // validate numeric arguments
  if (dt < 0)
  {
    std::cerr << "Time step must be a positive number!" << std::endl;
    return 7;
  }
  if (nskip < 0)
  {
    std::cerr << "Number of frames to skip must be at least 0!" << std::endl;
    return 8;
  }
  if (cutoff <= 0.0f)
  {
    std::cerr << "Cutoff distance must be greater than 0!" << std::endl;
    return 9;
  }
  if (mol_cutoff <= 0.0f)
  {
    std::cerr << "Molecule cutoff distance must be greater than 0!" << std::endl;
    return 10;
  }
  if (t_begin < 0.0f)
  {
    std::cerr << "Start time must be at least 0!" << std::endl;
    return 11;
  }
  if (t_end < t_begin && t_end != -1.0f)
  {
    std::cerr << "End time must be greater than start time!" << std::endl;
    return 12;
  }

  // create output directory if a prefix was given
  if (!out_prefix.empty())
  {
    std::error_code ec;
    std::filesystem::create_directories(std::filesystem::path(out_prefix), ec);
    if (ec && !std::filesystem::is_directory(out_prefix))
    {
      std::cerr << "Could not create output directory: " << out_prefix << std::endl;
      return 5;
    }
    if (!ec)
      std::cout << "WARNING: output directory already exists — files may be overwritten.\n";
  }

  cmdata::CMData cmdata(
    top_path, traj_path, cutoff, mol_cutoff, nskip, dt,
    mode, bkbn_H, weights_path, nopbc, t_begin, t_end, h5
  );
  cmdata.run();
  cmdata.process_data();
  cmdata.write_output(out_prefix);

  return 0;
}
