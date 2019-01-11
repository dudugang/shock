#include <Algebra.h>
#include <Eigen/Dense>
#include <gtest/gtest.h>
using namespace Eigen;

TEST(return_postive_test, ShouldReturnSameMatrixWhenPositive) {
    Matrix3d a;
    a << 1, 2, 3,
         4, 5, 6,
         7, 8, 9;
    Matrix3d b;
    b.setIdentity();

    ASSERT_EQ(a, Algebra::return_positive(a));
    ASSERT_EQ(b, Algebra::return_positive(b));
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
