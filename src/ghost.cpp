#include <iostream>
#include <ghost.h>
using std::vector;

// Constructor
Ghost::Ghost(double x, vector<double> q, int cell_id, vector<Cell*> neighbors,
    Face *left_face, Face *right_face) : Cell(x, q, cell_id, neighbors,
    left_face, right_face) {
    this->type = "ghost";
}


// Update ghost cells according to boundary conditions. As of now, reflective
// BCs are used.
void Ghost::update() {
    q[0] =  neighbors[0]->q[0];
    q[1] = -neighbors[0]->q[0];
    q[2] =  neighbors[0]->q[0];
}
