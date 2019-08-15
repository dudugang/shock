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


// Update ghost cell with an inflow boundary condition.
// TODO: Make this not hardcoded
void Ghost::update_inflow(unordered_map<int, Cell*> &cells) {

    double rho = 1.225;
    double u = 0;
    double v = 0;
    double p = 1e6;
    double gamma = 1.4;

    q[0] = rho;
    q[1] = rho * u;
    q[2] = rho * v;
    q[3] = p/(gamma-1) + (1/2)*rho*(u*u + v*v);

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
    q[3] =  neighbor->q[3];

    // Keep tangential momentum the same, but flip the sign of normal momentum.
    // This is actually somewhat nontrivial and I do not gaurantee that this is
    // correct.
    // TODO: Test this very well
    double &sintheta = faces[0]->sintheta;
    double &costheta = faces[0]->costheta;
    double &q1 = neighbor->q[1];
    double &q2 = neighbor->q[2];
    q[1] = 2*q2*costheta*sintheta + q1*(sintheta*sintheta - costheta*costheta);
    q[2] = q[1]*costheta/sintheta + q1*costheta/sintheta - q2;

}
