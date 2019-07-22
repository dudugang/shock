#include <cmath>
#include <iostream>
#include <vector>
#include <flowfield.h>
using std::cout;
using std::endl;
using std::vector;

Flowfield::Flowfield(Inputs inputs) {
    // Store inputs
    this->inputs = inputs;

    // Create initial conditions
    for (int i; i < inputs.n_cells; i++) {

        // Choose right or left of shock tube
        vector<double> q;
        q.reserve(3);
        // Shock
        if (i < inputs.n_cells/2) {
            q = inputs.q_left;
        } else {
            q = inputs.q_right;
        }

        // Find cell center
        double x = i*inputs.dx + inputs.dx/2;

        // Vector of cell neighbor IDs
        vector<int> neighbors{i-1, i+1};

        // Create cell and add to map of cells
        Cell current_cell = Cell(x, q, i, neighbors);
        cells[i] = &current_cell;

    }
    // Add ghost cells
    Ghost left_ghost  = Ghost(-inputs.dx/2, inputs.q_left, -1, vector<int>{0});
    Ghost right_ghost = Ghost(inputs.n_cells*inputs.dx + inputs.dx/2,
        inputs.q_right, inputs.n_cells, vector<int>{inputs.n_cells-1});

    // Create faces of cells and assign neighbor information
    for (int i = 0; i < inputs.n_cells + 1; i++) {
        Cell *left_cell = cells[i-1];
        Cell *right_cell = cells[i];
        Face current_face = Face(i, left_cell, right_cell);
        faces[i] = &current_face;
    }

}


// Calculate fluxes through every face in domain
void Flowfield::calculate_flux() {

    // For every pair of face IDs and faces, calculate the flux through each
    // of the faces
    for (auto &pair : faces) {

        // Convenient reference to the current face
        Face &face = *(pair.second);

        // Get face flux
        face.flux = Flux::steger_warming(face.left_cell->q, face.right_cell->q,
            inputs.gamma);
    }

}
/*

    def apply_time_integrator(self):
        for cell in self.cells.values():
            if cell.type == 'flow':
                cell.u = cell.u - (self.inputss.dt/self.inputss.dx)*(cell.flux_right - cell.flux_left)

    def update_ghosts(self):
        for cell in self.cells.values():
            if cell.type == 'ghost':
                cell.update(self.cells[cell.neighbors[0]])
*/
