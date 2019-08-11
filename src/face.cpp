#include <face.h>


// Constructor
Face::Face(Point point1, Point point2) {

    // Make sure that the points are oriented according to the sign conventions
    // illustrated in the header file
    if (point1.x > point2.x or point1.y > point2.y) {
        this->point1 = point1;
        this->point2 = point2;
    } else {
        this->point1 = point2;
        this->point2 = point1;
    }

    flux.resize(4);
    q_left.resize(4);
    q_right.resize(4);
    neighbors.reserve(2);

    // Find centroid of face
    center = Geometry::find_midpoint(point1, point2);

    // Calculate angle of normal vector, using sign conventions illustrated in
    // the header file for this class
    double dx = point1.x - point2.x;
    double dy = point1.y - point2.y;
    // Area is just the length of the line segment between the two points, since
    // this code is 2D (for now).
    area = std::sqrt(dx*dx + dy*dy);
    sintheta = dy/area;
    costheta = dx/area;
    theta = std::atan2(sintheta, costheta);
}


// Check if face nodes contain certain points
bool Face::contains(Point p1, Point p2) {

    // Output
    bool result;

    // Check if object contains the points in any order
    if ((point1 == p1) and (point2 == p2)) {
        result = true;
    } else if ((point1 == p2) and (point2 == p1)) {
        result = true;
    } else {
        result = false;
    }

    return result;

}
