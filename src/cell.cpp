#include <cell.h>
#include <iostream>
using std::vector;

// Constructor
Cell::Cell(vector<Point> vertices, vector<double> q, int volume_id,
    vector<Volume*> neighbors, Face *left_face, Face *right_face)
    : Volume(vertices, q, volume_id, neighbors, left_face, right_face) {
    this->type = "flow";
}

void Cell::update() {}
