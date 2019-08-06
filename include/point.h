#pragma once


// Class for a point in space. This can represent things like vertices in the
// mesh or cell centers.
class Point {

    public:
        Point(double, double);
        Point();
        double x;
        double y;

};
