#include <geometry.h>
#include <point.h>
#include <gtest/gtest.h>
using std::abs;


TEST(find_quad_centroid_test, ShouldGiveZeroForSquareAtOrigin) {

    // Create test data
    Point a(-1, -1);
    Point b(1, 1);
    Point c(1, -1);
    Point d(-1, 1);
    Point centroid = Geometry::find_quad_centroid(a, b, c, d);
    Point centroid_correct(0, 0);

    // Test
    double error = 1e-8;
    ASSERT_NEAR(centroid_correct.x, centroid.x, error);
    ASSERT_NEAR(centroid_correct.y, centroid.y, error);

}
