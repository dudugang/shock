#include <cell.h>
#include <iostream>
using std::vector;

// Constructor
Cell::Cell(vector<Point> vertices, vector<double> q, vector<Face*> faces,
    int volume_id) : Volume(vertices, q, faces, volume_id) {
    this->type = "flow";
}

void Cell::update() {}
