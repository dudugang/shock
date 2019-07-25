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
    for (int i = 0; i < inputs.n_cells; i++) {

        // Choose right or left of shock tube
        vector<double> q;
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
        Cell *current_cell = new Cell(x, q, i, neighbors, nullptr, nullptr);
        cells[i] = current_cell;

    }

    // Add ghost cells
    Ghost *left_ghost  = new Ghost(-inputs.dx/2, inputs.q_left, -1,
        vector<int>{0}, nullptr, nullptr);
    cells[-1] = left_ghost;
    Ghost *right_ghost = new Ghost(inputs.n_cells*inputs.dx + inputs.dx/2,
        inputs.q_right, inputs.n_cells, vector<int>{inputs.n_cells-1},
        nullptr, nullptr);
    cells[inputs.n_cells] = right_ghost;

    // Create faces of cells and assign neighbor information
    for (int i = 0; i < inputs.n_cells + 1; i++) {
        Cell *left_cell = cells[i-1];
        Cell *right_cell = cells[i];
        Face *current_face = new Face(i, left_cell, right_cell);
        // Inform neighbor cells of their new face
        left_cell->right_face = current_face;
        right_cell->left_face = current_face;
        faces[i] = current_face;
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


// Integrate semi-discrete equations in time
void Flowfield::apply_time_integrator() {

    // For every pair of cells and cell IDs, apply the time integrator
    for (auto &pair : cells) {

        // Convenient reference to the current cell
        Cell &cell = *(pair.second);

        // Only apply to flowfield cells, not ghosts
        if (cell.type == "flow") {
            // Apply to every equation
            for (int i = 0; i < inputs.n_equations; i++) {
                cell.q[i] = cell.q[i] - (inputs.dt/inputs.dx)
                    * (cell.right_face->flux[i] - cell.left_face->flux[i]);
            }
            // Debug output
            cout << cell.q[0] << " " << cell.q[1] << " " << cell.q[2] << endl;
        }

    }

}
/*
    def update_ghosts(self):
        for cell in self.cells.values():
            if cell.type == 'ghost':
                cell.update(self.cells[cell.neighbors[0]])
*/
