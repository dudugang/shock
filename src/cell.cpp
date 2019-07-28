#include <iostream>
#include <cell.h>
using std::vector;

// Constructor
Cell::Cell(double x, vector<double> q, int volume_id, vector<Volume*> neighbors,
    Face *left_face, Face *right_face) : Volume(x, q, volume_id, neighbors,
    left_face, right_face) {
    this->type = "flow";
}

// TODO: Add update, so that polymorphism could be used for both the time
// integrator and the ghost cell update
// TODO: Maybe don't do this, since flow cells need to be updated before ghosts
void Cell::update() {}
