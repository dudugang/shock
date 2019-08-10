#pragma once
#include <iostream>
#include <cmath>
#include <vector>
#include <point.h>
#include <volume.h>
using std::vector;


// Faces are 1D lines for 2D volumes and look something like this:
//
//                    O   <--- Point 1
//                   /                                y
//                  /
//                 /                                  |
//                /--\                                |
//               /    \--->   Normal vector           *----  x
//              /
//             /
//            O   <--- Point 2
//
// The normal vector always points to the right, unless the face is perfectly
// horizontal, in which case it points straight up. Point 1 is always
// counterclockwise of the normal vector, and point 2 is always clockwise of the
// normal vector. The angle theta is the angle between the x axis and the normal
// vector, with a positive theta denoting a normal vector above the x axis (the
// picture above shows a face with negative theta).
class Face {
    public:
        Face(Point, Point);
        vector<double> q_left;
        vector<double> q_right;
        Point point1;
        Point point2;
        double theta;
        double sintheta;
        double costheta;
        double ds;
        vector<double> flux;
};
