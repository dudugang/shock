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
        vector<Volume*> neighbors{id_to_volume[i-1], id_to_volume[i+1]};

        // Create cell and add to map/set of cells
        Cell *current_cell = new Cell(x, q, i, neighbors, nullptr, nullptr);
        id_to_volume[i] = current_cell;
        cells.insert(current_cell);
        volumes.insert(current_cell);

    }

    // Add ghost cells
    Ghost *left_ghost  = new Ghost(-inputs.dx/2, inputs.q_left, -1,
        vector<Volume*>{id_to_volume[0]}, nullptr, nullptr);
    id_to_volume[-1] = left_ghost;
    ghosts.insert(left_ghost);
    volumes.insert(left_ghost);
    Ghost *right_ghost = new Ghost(inputs.n_cells*inputs.dx + inputs.dx/2,
        inputs.q_right, inputs.n_cells, vector<Volume*>{id_to_volume[inputs.n_cells-1]},
        nullptr, nullptr);
    id_to_volume[inputs.n_cells] = right_ghost;
    ghosts.insert(right_ghost);
    volumes.insert(right_ghost);

    // Create faces of cells and assign neighbor information
    for (int i = 0; i < inputs.n_cells + 1; i++) {
        Volume *left_cell = id_to_volume[i-1];
        Volume *right_cell = id_to_volume[i];
        Face *current_face = new Face(i, left_cell, right_cell);
        // Inform neighbor cells of their new face
        left_cell->right_face = current_face;
        right_cell->left_face = current_face;
        id_to_face[i] = current_face;
        faces.insert(current_face);
    }

    // Start from time 0
    time = 0;

}


// Calculate fluxes through every face in domain
void Flowfield::calculate_flux(Flux &flux) {

    // For every pair of face IDs and faces, calculate the flux through each
    // of the faces
    for (auto &face : faces) {

        // Get face flux
        face->flux = flux.steger_warming(face->left_volume->q,
            face->right_volume->q, inputs.gamma);

    }

}


// Integrate semi-discrete equations in time
void Flowfield::apply_time_integrator() {

    // For every pair of cells and cell IDs, apply the time integrator
    for (auto &cell : cells) {

        // Apply to every equation
        for (int i = 0; i < inputs.n_equations; i++) {
            cell->q[i] = cell->q[i] - (inputs.dt/inputs.dx)
                * (cell->right_face->flux[i] - cell->left_face->flux[i]);
        }
        // Update time
        time = time + inputs.dt;

    }

    // Update ghost cells
    for (auto &ghost : ghosts) {
        ghost->update();
    }

}
