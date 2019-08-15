#pragma once
#include <algorithm>
#include <cmath>
#include <map>
#include <vector>
#include <point.h>
using std::atan2;
using std::map;
using std::vector;


// This class defines static functions for calculating geometric quantities,
// such as centroids, lengths, areas, or volumes.
class Geometry {

    public:
        static Point find_triangle_centroid(Point, Point, Point);
        static Point find_quad_centroid(Point, Point, Point, Point);
        static double find_quad_area(Point, Point, Point, Point);
        static double find_triangle_area(Point, Point, Point);
        static Point find_midpoint(Point, Point);
        static double angle_between(double, double);
        static constexpr double pi = 3.141592653589793238462643;

};
