#pragma once
#include <algorithm>
#include <cmath>
#include <iostream>
#include <unordered_map>
#include <vector>
#include <geometry.h>
#include <point.h>
using std::unordered_map;
using std::vector;


// Forward declare this since Face, Cell, and Ghost all depend on each other
class Volume;
class Cell;
class Ghost;

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
        void find_normal_vector(unordered_map<int, Cell*>&,
            unordered_map<int, Ghost*>&, int);
        bool contains(Point, Point);
        vector<double> q_left;
        vector<double> q_right;
        vector<int> neighbors;
        Point point1;
        Point point2;
        Point center;
        double theta;
        double sintheta;
        double costheta;
        double area;
        vector<double> flux;
};
