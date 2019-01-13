#include <Algebra.h>
#include <Eigen/Dense>
#include <gtest/gtest.h>
using namespace Eigen;

TEST(return_postive_test, ShouldReturnSameMatrixWhenPositive) {
    // Create test data
    Matrix3d a;
    a << 1, 2, 3,
         4, 5, 6,
         7, 8, 9;
    Matrix3d b;
    b.setIdentity();

    // Test
    ASSERT_EQ(a, Algebra::return_positive(a));
    ASSERT_EQ(b, Algebra::return_positive(b));
}

TEST(return_postive_test, ShouldReturnZeroMatrixWhenNegative) {
    // Create test data
    Matrix3d a;
    a << -1, -2, -3,
         -4, -5, -6,
         -7, -8, -9;
    Matrix3d b;
    b.setIdentity();
    b = -1*b;
    Matrix3d zero;
    zero.setZero();

    // Test
    ASSERT_EQ(zero, Algebra::return_positive(a));
    ASSERT_EQ(zero, Algebra::return_positive(b));
}

TEST(return_negative_test, ShouldReturnSameMatrixWhenNegative) {
    // Create test data
    Matrix3d a;
    a << -1, -2, -3,
         -4, -5, -6,
         -7, -8, -9;
    Matrix3d b;
    b.setIdentity();
    b = -1*b;

    // Test
    ASSERT_EQ(a, Algebra::return_negative(a));
    ASSERT_EQ(b, Algebra::return_negative(b));
}

TEST(return_negative_test, ShouldReturnZeroMatrixWhenPositive) {
    // Create test data
    Matrix3d a;
    a << 1, 2, 3,
         4, 5, 6,
         7, 8, 9;
    Matrix3d b;
    b.setIdentity();
    Matrix3d zero;
    zero.setZero();

    // Test
    ASSERT_EQ(zero, Algebra::return_negative(a));
    ASSERT_EQ(zero, Algebra::return_negative(b));
}
