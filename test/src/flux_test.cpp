#include <vector>
#include <flux.h>
#include <gtest/gtest.h>
using std::vector;


TEST(steger_warming_test, ShouldGiveZeroForZeroVelocityAndPressure) {

    // Create test data
    Flux flux;
    double rho = 1.225;
    double u = 0;
    double v = 0;
    double p = 0;
    double gamma = 1.4;
    vector<double> q_left  = {rho, rho*u, rho*v, p/(gamma-1) + (1/2)*rho*(u*u + v*v)};
    vector<double> q_right = {rho, rho*u, rho*v, p/(gamma-1) + (1/2)*rho*(u*u + v*v)};
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
    double rho = 1.225;
    double u = 0;
    double v = 0;
    double p = 1e6;
    double gamma = 1.4;
    vector<double> q_left  = {rho, rho*u, rho*v, p/(gamma-1) + (1/2)*rho*(u*u + v*v)};
    vector<double> q_right = {rho, rho*u, rho*v, p/(gamma-1) + (1/2)*rho*(u*u + v*v)};
    vector<double> face_flux = flux.steger_warming(q_left, q_right, gamma);
    vector<double> face_flux_correct = {0, p, 0, 0};

    // Test
    double error = 1e-8;
    ASSERT_NEAR(face_flux_correct[0], face_flux[0], error);
    ASSERT_NEAR(face_flux_correct[1], face_flux[1], error);
    ASSERT_NEAR(face_flux_correct[2], face_flux[2], error);
    ASSERT_NEAR(face_flux_correct[3], face_flux[3], error);

}
