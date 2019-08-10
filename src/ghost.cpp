#include <ghost.h>


// Constructor
Ghost::Ghost(vector<Point> vertices, vector<double> q, vector<Face*> faces,
    int ghost_id) : Volume(vertices, q, faces, ghost_id) {
    this->type = "ghost";
    cout << "Creating ghost with ID " << volume_id << endl;
}


// Update ghost cells according to boundary conditions. As of now, reflective
// BCs are used.
// TODO: Make this work in 2D in arbitrary directions. Must use some info from
// normal vector angle for the momentum equations.
void Ghost::update() {
    /*
    q[0] =  neighbors[0]->q[0];
    q[1] = -neighbors[0]->q[1];
    q[2] = -neighbors[0]->q[1];
    q[3] =  neighbors[0]->q[2];
    TODO: fix this
    */
}
