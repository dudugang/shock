#include <cmath>
#include <geometry.h>
#include <point.h>
#include <gtest/gtest.h>
using std::abs;
using std::sqrt;


TEST(find_triangle_centroid_test, ShouldGiveZeroForTriangleAtOrigin) {

    // Create test data
    // This is an equilateral triangle, side lengths equal to 1, with centroid
    // lining up with the origin because of the 2/3-rule for triangle centroids
    Point a(0, sqrt(3)/3);
    Point b(.5, -sqrt(3)/6);
    Point c(-.5, -sqrt(3)/6);
    Point centroid = Geometry::find_triangle_centroid(a, b, c);
    Point centroid_correct(0, 0);

    // Test
    double error = 1e-8;
    ASSERT_NEAR(centroid_correct.x, centroid.x, error);
    ASSERT_NEAR(centroid_correct.y, centroid.y, error);

}


TEST(find_quad_centroid_test, ShouldGiveZeroForSquareAtOrigin) {

    // Create test data
    // This is a square, with its center at the origin
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


TEST(find_triangle_area_test, ShouldMatchAreaOfRightTriangle) {

    // Create test data
    // This is a right triangle, the area of which can be found with (1/2)*b*h
    Point a(-1, -1);
    Point b(2, -1);
    Point c(2, 3);
    double area = Geometry::find_triangle_area(a, b, c);
    double area_correct = 6;

    // Test
    double error = 1e-8;
    ASSERT_NEAR(area_correct, area, error);

}


TEST(find_quad_area_test, ShouldMatchAreaOfRectangle) {

    // Create test data
    // This is a rectangle, the area of which can be found with b*h
    Point a(-1, -1);
    Point b(3, -1);
    Point c(3, 2);
    Point d(-1, 2);
    double area = Geometry::find_quad_area(a, b, c, d);
    double area_correct = 12;

    // Test
    double error = 1e-8;
    ASSERT_NEAR(area_correct, area, error);

}


TEST(find_midpoint_test, ShouldGiveZeroForLineAroundOrigin) {

    // Create test data
    // This is a line, with endpoints equally far from the origin on either side
    Point a(2, 1);
    Point b(-2, -1);
    Point midpoint = Geometry::find_midpoint(a, b);
    Point midpoint_correct(0, 0);

    // Test
    double error = 1e-8;
    ASSERT_NEAR(midpoint_correct.x, midpoint.x, error);
    ASSERT_NEAR(midpoint_correct.y, midpoint.y, error);

}


TEST(angle_between_test, ShouldGiveAngleBetweenForEasyAngles) {

    // Create test data
    double alpha1 = 0;
    double beta1  = 0;
    double alpha2 = 0;
    double beta2  = .75*Geometry::pi;
    double alpha3 = (1/6)*Geometry::pi;
    double beta3  = (2/3)*Geometry::pi;
    double difference1 = Geometry::angle_between(alpha1, beta1);
    double difference2 = Geometry::angle_between(alpha2, beta2);
    double difference3 = Geometry::angle_between(alpha3, beta3);
    double difference_correct1 = 0;
    double difference_correct2 = .75*Geometry::pi;
    double difference_correct3 = (7/12)*Geometry::pi;

    // Test
    double error = 1e-8;
    ASSERT_NEAR(difference_correct1, difference1, error);
    ASSERT_NEAR(difference_correct2, difference2, error);
    ASSERT_NEAR(difference_correct3, difference3, error);

}
