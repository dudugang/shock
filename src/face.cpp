#include <face.h>
using std::cout;
using std::endl;


// Constructor
Face::Face(Point point1, Point point2) {
    this->point1 = point1;
    this->point2 = point2;
    flux.resize(3);
    q_left.resize(4);
    q_right.resize(4);

    // Calculate angle of normal vector, using sign conventions illustrated in
    // the header file for this class
    double dx = point1.x - point2.x;
    double dy = point1.y - point2.y;
    ds = std::sqrt(dx*dx + dy*dy);
    sintheta = dy/ds;
    costheta = dx/ds;
    theta = std::atan2(sintheta, costheta);
}
