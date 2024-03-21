#include "nbody.h"
#include <gtest/gtest.h>
#include "vec3.h"

class NBodyTest : public testing::Test{
  protected:
    void SetUp() override {
      std::vector<Body> particles;
      particles.emplace_back(6.66666667,
                             vec3(0.27626589, -1.85462808, 0.62390111),
                             vec3(-0.43778315, 2.171257, 1.15231025));
      particles.emplace_back(6.66666667,
                             vec3(1.14531129, 1.03719047, 1.88663893),
                             vec3(-1.81881234, -0.13804934, 0.53983961));
      particles.emplace_back(6.66666667,
                             vec3(-0.11169829, -0.36210134, 0.14867505),
                             vec3(-1.77528229, 1.31487654, -0.47344805));
      simulation = {particles, 1};
    }

    void TearDown() override{
    }

    NBody simulation;

};

TEST_F(NBodyTest, SuperMegaHardTest){
  ASSERT_EQ(simulation.gravitationalConstant, 1);
}

TEST_F(NBodyTest, AllPairsTestPosition){
  simulation.update(0.01);
  std::vector<Body> estimated = simulation.GetParticles();
  std::vector<vec3>  trueValues = {
      {0.27626589, -1.85462808, 0.62390111},
      {1.14531129, 1.03719047, 1.88663893},
      {-0.11169829, -0.36210134, 0.14867505}
  };
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      ASSERT_NEAR(estimated[i].GetPosition()[j], trueValues[i][j], 0.1);
    }
  }
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}