#include <vector>
#include <flux.h>
#include <state.h>
#include <gtest/gtest.h>
using std::vector;


TEST(steger_warming_test, ShouldGiveZeroForZeroVelocityAndPressure) {

    // Create test data
    Flux flux;
    double gamma = 1.4;
    double rho = 1.225;
    double u = 0;
    double v = 0;
    double p = 0;
    vector<double> physical = {rho, u, v, p};
    vector<double> q_left  = State::physical_to_conserved(physical, gamma);
    vector<double> q_right = q_left;
    vector<double> face_flux = flux.steger_warming(q_left, q_right, gamma);
    vector<double> face_flux_correct = {0, 0, 0, 0};

    // Test
    double error = 1e-8;
    ASSERT_NEAR(face_flux_correct[0], face_flux[0], error);
    ASSERT_NEAR(face_flux_correct[1], face_flux[1], error);
    ASSERT_NEAR(face_flux_correct[2], face_flux[2], error);
    ASSERT_NEAR(face_flux_correct[3], face_flux[3], error);

}


TEST(steger_warming_test, ShouldGivePressureForZeroVelocity) {

    // Create test data
    Flux flux;
    double gamma = 1.4;
    double rho = 1.225;
    double u = 0;
    double v = 0;
    double p = 1e6;
    vector<double> physical = {rho, u, v, p};
    vector<double> q_left  = State::physical_to_conserved(physical, gamma);
    vector<double> q_right = q_left;
    vector<double> face_flux = flux.steger_warming(q_left, q_right, gamma);
    vector<double> face_flux_correct = {0, p, 0, 0};

    // Test
    double error = 1e-8;
    ASSERT_NEAR(face_flux_correct[0], face_flux[0], error);
    ASSERT_NEAR(face_flux_correct[1], face_flux[1], error);
    ASSERT_NEAR(face_flux_correct[2], face_flux[2], error);
    ASSERT_NEAR(face_flux_correct[3], face_flux[3], error);

}


TEST(steger_warming_test, ShouldBeUpwindedForSupersonicFlow) {

    // Create test data
    Flux flux;
    double gamma = 1.4;
    double rho = 1;
    double u = 1000;
    double v = 0;
    double p = 1e5;
    vector<double> physical = {rho, u, v, p};
    vector<double> q_left  = State::physical_to_conserved(physical, gamma);
    vector<double> q_right = q_left;
    flux.steger_warming(q_left, q_right, gamma);

    // Test
    for (int i = 0; i < q_left.size(); i++) {
        ASSERT_TRUE(flux.eig_plus[i] >= 0);
        ASSERT_TRUE(flux.eig_minus[i] == 0);
        ASSERT_TRUE(flux.flux_plus[i] >= 0);
        ASSERT_TRUE(flux.flux_minus[i] == 0);
    }

}
