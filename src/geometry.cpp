#include <geometry.h>


// Find the centroid of a triangle from its vertices
Point Geometry::find_triangle_centroid(Point a, Point b, Point c) {

    // Find centroid. Averaging the vertices works for triangles and quads.
    return (a + b + c) / 3.0;

}


// Find the centroid of a quadrilateral from its vertices
Point Geometry::find_quad_centroid(Point a, Point b, Point c, Point d) {

    // Find centroid. Averaging the vertices works for triangles and quads.
    return (a + b + c + d) / 4.0;

}


// Find the area of a triangle from its vertices
double Geometry::find_triangle_area(Point a, Point b, Point c) {

    // Find area. This is essentially a cross product of two of the triangle's
    // sides divided by 2, since the magnitude of the cross product of two
    // vectors equals the area of the parallelogram between the vectors, and the
    // triangle is just half the parallelogram.
    return (a.x*(b.y - c.y) + b.x*(c.y - a.y) + c.x*(a.y - b.y)) / 2;

}


// Find the area of a quadrilateral from its vertices
double Geometry::find_quad_area(Point a, Point b, Point c, Point d) {

    // Find area. This amounts to finding the areas of the two triangles that
    // form the quadrilateral and adding them together. However, there is an
    // extra complication since the vertices may not be ordered, so we need to
    // order them first.

    // Sort cells in clockwise order, assuming convex cells
    // TODO: Make sure code checks to make sure no cells are concave
    // Find any point inside polygon, in this case the centroid
    Point centroid = Geometry::find_quad_centroid(a, b, c, d);

    // Find angles between centroid and all points
    vector<double> angles;
    angles.resize(4);
    angles[0] = atan2(a.y - centroid.y, a.x - centroid.x);
    angles[1] = atan2(b.y - centroid.y, b.x - centroid.x);
    angles[2] = atan2(c.y - centroid.y, c.x - centroid.x);
    angles[3] = atan2(d.y - centroid.y, d.x - centroid.x);

    // Map angles to points
    map<double, Point> angles_to_points;
    angles_to_points[angles[0]] = a;
    angles_to_points[angles[1]] = b;
    angles_to_points[angles[2]] = c;
    angles_to_points[angles[3]] = d;

    // Sort in increasing order
    std::sort(angles.begin(), angles.end());

    // Find areas of two triangles that make up the quadrilateral using the now
    // sorted vertices
    double area1 = Geometry::find_triangle_area(angles_to_points[angles[0]],
        angles_to_points[angles[1]], angles_to_points[angles[2]]);
    double area2 = Geometry::find_triangle_area(angles_to_points[angles[0]],
        angles_to_points[angles[2]], angles_to_points[angles[3]]);

    // Add two triangle areas to get quadrilateral area
    return area1 + area2;

}


// Find the midpoint of a line from its end points
Point Geometry::find_midpoint(Point a, Point b) {

    // Find midpoint. Averaging end points works for straight lines.
    return (a + b) / 2.0;

}


// Find the angle between two vectors defined by their angles from a common
// axis.
double Geometry::angle_between(double alpha, double beta) {

    // Get difference in angles
    double difference = std::abs(alpha - beta);

    // Difference can not be greater than pi radians, since after pi radians of
    // rotation the two vectors start getting close to each other, not further
    // away.
    if (difference > pi) {
        difference = difference - 2*(difference - pi);
    }

}
