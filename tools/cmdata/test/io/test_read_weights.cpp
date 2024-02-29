#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "io.hpp"

class ReadWeightsFileTest : public ::testing::Test {
protected:
};

TEST(ReadWeightsFileTest, CorrectFileContents)
{
  std::vector<double> weights = cmdata::io::read_weights_file("resources/weights_good.txt");
  ASSERT_THAT(weights, testing::ElementsAre(
    5.20805e-02, 6.26943e-02, 5.89982e-02, 5.90832e-02, 6.20511e-02,
    1.49853e-02, 3.71899e-02, 7.67860e-02, 1.94643e-02, 5.78728e-02,
    3.24808e-02, 6.95297e-02, 7.23627e-02, 6.87851e-02, 3.19701e-02,
    3.56076e-02, 3.26592e-02, 6.14780e-02, 7.14015e-02, 2.25196e-02
  ));
}

TEST(ReadWeightsFileTest, FileNotFound)
{
  std::string filename = "resources/weights_nofile.txt";
  ASSERT_THROW(cmdata::io::read_weights_file(filename), std::runtime_error);
}

TEST(ReadWeightsFileTest, EmptyFile)
{
  std::string filename = "resources/weights_empty.txt";
  ASSERT_THROW(cmdata::io::read_weights_file(filename), std::runtime_error);
}

TEST(ReadWeightsFileTest, NegativeValues)
{
  std::string filename = "resources/weights_bad.txt";
  ASSERT_THROW(cmdata::io::read_weights_file(filename), std::runtime_error);
}
