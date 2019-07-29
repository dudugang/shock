#include <iostream>
#include <ghost.h>
using std::vector;

// Constructor
Ghost::Ghost(vector<double> center, vector<double> q, int volume_id,
    vector<Volume*> neighbors, Face *left_face, Face *right_face)
    : Volume(center, q, volume_id, neighbors, left_face, right_face) {
    this->type = "ghost";
}


// Update ghost cells according to boundary conditions. As of now, reflective
// BCs are used.
// TODO: Make this work in 2D in arbitrary directions. Must use some info from
// normal vector angle for the momentum equations.
void Ghost::update() {
    q[0] =  neighbors[0]->q[0];
    q[1] = -neighbors[0]->q[1];
    q[2] = -neighbors[0]->q[1];
    q[3] =  neighbors[0]->q[2];
}
