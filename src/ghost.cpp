#include <ghost.h>
#include <bc.h>
#include <cell.h>
#include <face.h>


// Constructor
Ghost::Ghost(vector<Point> vertices, vector<double> q, vector<Face*> faces,
    int id) : Volume(vertices, q, faces, id) {
    this->type = "ghost";
    cout << "Creating ghost with ID " << id << endl;
}


// Update ghost cell according to its boundary condition type
void Ghost::update(Inputs &inputs, unordered_map<int, Cell*> &cells) {

    BC boundary_condition = inputs.bc[type];
    switch(boundary_condition) {
        case BC::inflow:
            update_inflow(cells);
            break;
        case BC::wall:
            update_wall(cells);
            break;
    }

}


void Ghost::update_inflow(unordered_map<int, Cell*> &cells) {
}


// Update ghost cell with an inviscid wall boundary condition.
// TODO: Make this work in 2D in arbitrary directions. Must use some info from
// normal vector angle for the momentum equations.
void Ghost::update_wall(unordered_map<int, Cell*> &cells) {

    // Find neighbor. A ghost only has one face, so faces[0] retrieves this, and
    // ghosts always have larger cell ID's than flowfield cells, so neighbors[0]
    // retreives the neighboring cell ID since neighbors is sorted from small
    // cell ID to large cell ID.
    Cell* neighbor = cells[faces[0]->neighbors[0]];

    // Keep mass and energy the same as its neighbor
    q[0] =  neighbor->q[0];
    q[3] =  neighbor->q[2];

    q[1] = -neighbor->q[1];
    q[2] = -neighbor->q[1];

}
