#include <point.h>


// Constructor
Point::Point(double x, double y) {

    this->x = x;
    this->y = y;

}

// Default constructor: do nothing
Point::Point () {}

// Check if two points equal
bool Point::operator== (const Point &p) {
    return (x == p.x and y == p.y);
}
