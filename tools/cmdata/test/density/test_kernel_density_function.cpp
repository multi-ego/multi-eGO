#include <gtest/gtest.h>
#include <vector>
#include <cmath>

#include <iostream>

#include "density.hpp"
#include "indexing.hpp"

class KernelDensityEstimatorTest : public ::testing::Test {
protected:
    void SetUp() override {}

    void TearDown() override {}
};

TEST_F(KernelDensityEstimatorTest, TestFunctionality)
{
  std::vector<double> bins = {0.01,0.02,0.03,0.04,0.05};
  std::vector<double> x = {0.0, 0.0, 0.0, 0.0, 0.0};
  double mu = 0.03;
  double norm = 1.0;
  
  cmdata::density::kernel_density_estimator(x.begin(), bins, mu, norm);  
  
  ASSERT_NEAR(x[0], 0.0, 0.0001);
  ASSERT_NEAR(x[1], 25.453, 0.0001);
  ASSERT_NEAR(x[2], 46.7075, 0.0001);
  ASSERT_NEAR(x[3], 25.453, 0.0001);
  ASSERT_NEAR(x[4], 0.0, 0.0001);
  
  bins = {-0.04, -0.03, -0.02, -0.01, 0.0, 0.01};
  x = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  mu = -0.02;
  norm = 1.0;
  
  cmdata::density::kernel_density_estimator(x.begin(), bins, mu, norm);  

  ASSERT_NEAR(x[0], 0.0, 0.0001);
  ASSERT_NEAR(x[1], 50.9061, 0.0001);
  ASSERT_NEAR(x[2], 93.4149, 0.0001);
  ASSERT_NEAR(x[3], 50.9061, 0.0001);
  ASSERT_NEAR(x[4], 0.0, 0.0001);
  ASSERT_NEAR(x[5], 0.0, 0.0001);
}

TEST_F(KernelDensityEstimatorTest, NullTest)
{
  std::vector<double> bins = {1.0, 2.0, 3.0, 4.0, 5.0};
  std::vector<double> x = {0.0, 0.0, 0.0, 0.0, 0.0};
  double mu = 9.0;
  double norm = 1.0;
  
  cmdata::density::kernel_density_estimator(x.begin(), bins, mu, norm);  
  
  ASSERT_EQ(x[0], 0.0);
  ASSERT_EQ(x[1], 0.0);
  ASSERT_EQ(x[2], 0.0);
  ASSERT_EQ(x[3], 0.0);
  ASSERT_EQ(x[4], 0.0);  
}

TEST_F(KernelDensityEstimatorTest, BorderCaseTest)
{
  std::vector<double> bins = {0.0, 0.005, 0.01, 0.015, 0.02};
  std::vector<double> x = {0.0, 0.0, 0.0, 0.0, 0.0};
  double mu = 0.0;
  double norm = 10.0;
  
  cmdata::density::kernel_density_estimator(x.begin(), bins, mu, norm);

  ASSERT_NEAR(x[0], 934.149, 0.001);
  ASSERT_NEAR(x[1], 807.203, 0.001);
  ASSERT_NEAR(x[2], 509.061, 0.001);
  ASSERT_NEAR(x[3], 204.531, 0.001);
  ASSERT_NEAR(x[4], 0.0, 0.001);
}