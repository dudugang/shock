#pragma once


// Class for a point in space. This can represent things like vertices in the
// mesh or cell centers. An ID is assigned to the point during construction.
class Point {

    public:
        Point(double, double, int);
        Point();
        double x;
        double y;
        int point_id;

};
