#pragma once


// Class for a point in space. This can represent things like vertices in the
// mesh or cell centers.
class Point {

    public:
        Point(double, double);
        Point();
        bool operator== (const Point &p);
        double x;
        double y;

};
