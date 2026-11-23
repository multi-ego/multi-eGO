#ifndef _CMDATA_WEIGHTS_HPP
#define _CMDATA_WEIGHTS_HPP

// Per-frame weight file parsing. Pure C++17, no GROMACS dependency, so it can
// be unit-tested on its own (see test/test_cmdata.cpp).
//
// Format: plain text, one non-negative finite weight per trajectory frame,
// one value per line, in frame order. Weights are linear (not log) and need
// not be normalised. Blank lines and lines whose first non-blank character is
// '#' or '@' (PLUMED / xvg headers) are ignored. Windows line endings are
// accepted. Any other line must contain exactly one number.

#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <istream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace cmdata::io
{

inline std::vector<float> parse_weights( std::istream &in, const std::string &source = "weights file" )
{
  std::vector<float> w;
  std::string line;
  std::size_t lineno = 0;

  auto fail = [&](const std::string &what)
  {
    throw std::runtime_error(source + ", line " + std::to_string(lineno) + ": " + what);
  };

  while (std::getline(in, line))
  {
    ++lineno;
    std::istringstream iss(line);
    std::string token;
    if (!(iss >> token)) continue;                   // blank or whitespace-only (incl. "\r")
    if (token[0] == '#' || token[0] == '@') continue; // comment / header

    std::string extra;
    if (iss >> extra)
      fail("expected one value per line, found more than one column (\"" + line +
           "\"); extract the weight column first");

    errno = 0;
    char *end = nullptr;
    double v = std::strtod(token.c_str(), &end);
    if (end == token.c_str() || *end != '\0')
      fail("cannot parse \"" + token + "\" as a number");
    if (errno == ERANGE || !std::isfinite(v) || !std::isfinite(static_cast<float>(v)))
      fail("weight \"" + token + "\" is not finite (or overflows float)");
    if (v < 0.0)
      fail("negative weight " + token);

    w.push_back(static_cast<float>(v));
  }

  if (w.empty())
    throw std::runtime_error(source + " contains no weights");

  return w;
}

inline std::vector<float> read_weights_file( const std::string &path )
{
  std::ifstream infile(path);
  if (!infile.good())
    throw std::runtime_error("Cannot open weights file " + path);
  return parse_weights(infile, "Weights file " + path);
}

} // namespace cmdata::io

#endif // _CMDATA_WEIGHTS_HPP
