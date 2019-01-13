#include <Flux.h>
#include <Eigen/Dense>
#include <gtest/gtest.h>
using namespace Eigen;
using std::abs;

TEST(calculate_f_right_test, ShouldGiveSameResultsAsMatlab) {
    // Create test data
    double gamma = 1.4;
    Vector3d a;
    Vector3d a_right;
    a << 1,
         1e-15,
         100;
    a_right << 1,
               1e-15,
               100;
    Vector3d a_f = Flux::calculate_f_right(a, a_right, gamma);
    Vector3d a_f_correct;
    a_f_correct << 8.8817841970012523e-16, // NOTE: This may or may not be right
                   40,
                   1.13686837721616e-13;
    Vector3d b;
    Vector3d b_right;
    b << 1.5,
         20,
         3700000;
    b_right << 1,
               10,
               2500000;
    Vector3d b_f = Flux::calculate_f_right(b, b_right, gamma);
    Vector3d b_f_correct;
    b_f_correct << 223.470064319141,
                   1248508.7353447,
                   747352019.685737;
    double tolerance = 1e-8;
    double error;

    // Test
    for (int i = 0; i < 3; i++) {
        error = abs(tolerance * a_f_correct(i));
        ASSERT_NEAR(a_f_correct(i), a_f(i), error) << "Failure occured at index i = " << i;
        error = abs(tolerance * b_f_correct(i));
        ASSERT_NEAR(b_f_correct(i), b_f(i), error) << "Failure occured at index i = " << i;
    }
}

TEST(calculate_f_left_test, ShouldGiveSameResultsAsMatlab) {
    // Create test data
    double gamma = 1.4;
    Vector3d a;
    Vector3d a_left;
    a << 1,
         1e-15,
         100;
    a_left << 1,
              1e-15,
              100;
    Vector3d a_f = Flux::calculate_f_left(a, a_left, gamma);
    Vector3d a_f_correct;
    a_f_correct << 8.8817841970012523e-16, // NOTE: This may or may not be right
                   40,
                   1.13686837721616e-13;
    Vector3d b;
    Vector3d b_left;
    b << 1,
         10,
         2500000;
    b_left << 1.5,
              20,
              3700000;
    Vector3d b_f = Flux::calculate_f_left(b, b_left, gamma);
    Vector3d b_f_correct;
    b_f_correct << 223.470064319141,
                   1248508.7353447,
                   747352019.685737;
    double tolerance = 1e-8;
    double error;

    // Test
    for (int i = 2; i < 3; i++) {
        error = abs(tolerance * a_f_correct(i));
        ASSERT_NEAR(a_f_correct(i), a_f(i), error) << "Failure occured at index i = " << i;
        error = abs(tolerance * b_f_correct(i));
        ASSERT_NEAR(b_f_correct(i), b_f(i), error) << "Failure occured at index i = " << i;
    }
}

TEST(calculate_right_eigenvectors_test, ShouldGiveSameResultsAsMatlab) {
    // Create test data
    double gamma = 1.4;
    Vector3d a;
    a << 1,
         1e-15,
         100;
    Matrix3d a_eigenvectors = Flux::calculate_right_eigenvectors(a, gamma);
    Matrix3d a_eigenvectors_correct;
    a_eigenvectors_correct << 2e+30, 0.00714285714285714, 0.00714285714285714,
                              2e+15, 0.0534522483824849, -0.0534522483824849,
                              1, 1, 1;
    Vector3d b;
    b << 1.5,
         20,
         3700000;
    Matrix3d b_eigenvectors = Flux::calculate_right_eigenvectors(b, gamma);
    Matrix3d b_eigenvectors_correct;
    b_eigenvectors_correct << 0.01125, 2.88270156418778e-07, 2.90898311798352e-07,
                              0.15, 0.000342641810039285, -0.000338008382784403,
                              1, 1, 1;
    double tolerance = 1e-8;
    double error;

    // Test
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            error = abs(tolerance * a_eigenvectors_correct(i,j));
            ASSERT_NEAR(a_eigenvectors_correct(i,j), a_eigenvectors(i,j), error);
            error = abs(tolerance * b_eigenvectors_correct(i,j));
            ASSERT_NEAR(b_eigenvectors_correct(i,j), b_eigenvectors(i,j), error);
        }
    }
}

TEST(calculate_left_eigenvectors_test, ShouldGiveSameResultsAsMatlab) {
    // Create test data
    double gamma = 1.4;
    Vector3d a;
    a << 1,
         1e-15,
         100;
    Matrix3d a_eigenvectors = Flux::calculate_left_eigenvectors(a, gamma);
    Matrix3d a_eigenvectors_correct;
    a_eigenvectors_correct << 5e-31, 3.57142857142857e-48, -3.57142857142857e-33,
                              -9.35414346693485e-15, 9.35414346693485, 0.5,
                              9.35414346693485e-15, -9.35414346693485, 0.5;
    Vector3d b;
    b << 1.5,
         20,
         3700000;
    Matrix3d b_eigenvectors = Flux::calculate_left_eigenvectors(b, gamma);
    Matrix3d b_eigenvectors_correct;
    b_eigenvectors_correct << 88.8866008041474, 0.000343212711225973, -2.5740953341948e-05,
                              -19632.7492120062, 1469.10764534752, 0.50228183294185,
                              19543.8626112021, -1469.10798856023, 0.497743908011492;
    double tolerance = 1e-8;
    double error;

    // Test
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            error = abs(tolerance * a_eigenvectors_correct(i,j));
            ASSERT_NEAR(a_eigenvectors_correct(i,j), a_eigenvectors(i,j), error);
            error = abs(tolerance * b_eigenvectors_correct(i,j));
            ASSERT_NEAR(b_eigenvectors_correct(i,j), b_eigenvectors(i,j), error);
        }
    }
}

TEST(calculate_eigenvalues_test, ShouldGiveSameResultsAsMatlab) {
    // Create test data
    double gamma = 1.4;
    Vector3d a;
    a << 1,
         1e-15,
         100;
    Matrix3d a_eigenvalues = Flux::calculate_eigenvalues(a, gamma);
    Matrix3d a_eigenvalues_correct;
    a_eigenvalues_correct << 1e-15, 0, 0,
                             0, 7.48331477354788, 0,
                             0, 0, -7.48331477354788;
    Vector3d b;
    b << 1.5,
         20,
         3700000;
    Matrix3d b_eigenvalues = Flux::calculate_eigenvalues(b, gamma);
    Matrix3d b_eigenvalues_correct;
    b_eigenvalues_correct << 13.3333333333333, 0, 0,
                             0, 1188.61353632986, 0,
                             0, 0, -1161.94686966319;
    double tolerance = 1e-8;
    double error;

    // Test
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            error = abs(tolerance * a_eigenvalues_correct(i,j));
            ASSERT_NEAR(a_eigenvalues_correct(i,j), a_eigenvalues(i,j), error);
            error = abs(tolerance * b_eigenvalues_correct(i,j));
            ASSERT_NEAR(b_eigenvalues_correct(i,j), b_eigenvalues(i,j), error);
        }
    }
}
