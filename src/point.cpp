#include <point.h>


// Constructor
Point::Point(double x, double y) {

    this->x = x;
    this->y = y;

}


// Default constructor: do nothing
Point::Point() {}


// Check if two points equal
bool operator== (const Point &p1, const Point &p2) {
    return (p1.x == p2.x and p1.y == p2.y);
}


// Add points
Point operator+ (const Point &p1, const Point &p2) {
    return Point(p1.x + p2.x, p1.y + p2.y);
}


// Subtract points
Point operator- (const Point &p1, const Point &p2) {
    return Point(p1.x - p2.x, p1.y - p2.y);
}


// Multiply point and double
Point operator* (const Point &p, const double &num) {
    return Point(p.x * num, p.y * num);
}


// Divide point and double
Point operator/ (const Point &p, const double &num) {
    return Point(p.x / num, p.y / num);
}
