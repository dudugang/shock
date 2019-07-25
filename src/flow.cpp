#include <iostream>
#include <flow.h>
using std::vector;

// Constructor
Flow::Flow(double x, vector<double> q, int cell_id, vector<Cell*> neighbors,
    Face *left_face, Face *right_face) : Cell(x, q, cell_id, neighbors,
    left_face, right_face) {
    this->type = "flow";
}

// TODO: Add update, so that polymorphism could be used for both the time
// integrator and the ghost cell update
void Flow::update() {}
