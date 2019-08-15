#include <vector>
#include <flux.h>
#include <gtest/gtest.h>
using std::vector;


TEST(steger_warming_test, ShouldGiveZeroForEqualLeftAndRightStates) {

    // Create test data
    Flux flux;
    vector<double> q_left = {1.225, 0, 0, 1e6};
    vector<double> q_right = {1.225, 0, 0, 1e6};
    double gamma = 1.4;
    vector<double> face_flux = flux.steger_warming(q_left, q_right, gamma);
    vector<double> face_flux_correct = {0, 0, 0, 0};

    // Test
    double error = 1e-8;
    //ASSERT_NEAR(face_flux_correct[0], face_flux[0], error);
    ASSERT_NEAR(face_flux_correct[1], face_flux[1], error);
    ASSERT_NEAR(face_flux_correct[2], face_flux[2], error);
    ASSERT_NEAR(face_flux_correct[3], face_flux[3], error);

}
