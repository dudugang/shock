#include <geometry.h>


// Find the centroid of a triangle from its vertices
Point Geometry::find_triangle_centroid(Point a, Point b, Point c) {

    // Find centroid
    return (a + b + c) / 3.0;

}


// Find the centroid of a quadrilateral from its vertices
Point Geometry::find_quad_centroid(Point a, Point b, Point c, Point d) {

    // Sort cells in clockwise order, assuming convex cells
    // TODO: Make sure code checks to make sure no cells are concave
    // Find any point inside polygon, in this case the average of all points
    Point average = (a + b + c + d) / 4;

    // Find angles between average and all points
    vector<double> angles;
    angles.resize(4);
    angles[0] = atan2(a.y - average.y, a.x - average.x);
    angles[1] = atan2(b.y - average.y, b.x - average.x);
    angles[2] = atan2(c.y - average.y, c.x - average.x);
    angles[3] = atan2(d.y - average.y, d.x - average.x);

    // Map angles to points
    map<double, Point> angles_to_points;
    angles_to_points[angles[0]] = a;
    angles_to_points[angles[1]] = b;
    angles_to_points[angles[2]] = c;
    angles_to_points[angles[3]] = d;

    // Sort in increasing order
    std::sort(angles.begin(), angles.end());

    // Find midpoints of opposite points
    Point midpoint1 = (angles_to_points[angles[0]] + angles_to_points[angles[2]]) / 2;
    Point midpoint2 = (angles_to_points[angles[1]] + angles_to_points[angles[3]]) / 2;

    // Find midpoint of midpoints
    return (midpoint1 + midpoint2) / 2;

}
