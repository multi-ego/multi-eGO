#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "io.hpp"

class ReadSymmetryIndicesTest : public ::testing::Test {
protected:
  std::string top_path = "resources/test_gp.tpr";
  gmx_mtop_t* mtop;
  int natoms;
  std::vector<int> natmol2;
  std::vector<int> start_index = {0};

  void SetUp() override
  {
    // read topology
    matrix boxtop;
    mtop = (gmx_mtop_t*)malloc(sizeof(gmx_mtop_t));
    TpxFileHeader header = readTpxHeader(top_path.c_str(), true);
    PbcType pbcType = read_tpx(top_path.c_str(), nullptr, boxtop, &natoms, nullptr, nullptr, mtop);
    natmol2.push_back(natoms);
  }

  void TearDown() override
  {
    free(mtop);
  }
};

TEST_F(ReadSymmetryIndicesTest, ThrowsExceptionForNonexistentFile)
{
  std::vector<std::vector<std::vector<int>>> eq_list;
  // Expect that the function throws a runtime_error when file does not exist
  EXPECT_THROW(
      cmdata::io::read_symmetry_indices("nonexistent_file.txt", mtop, eq_list, natmol2, start_index),
      std::runtime_error
  );
}