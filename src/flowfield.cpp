#include <flowfield.h>
#include <cmath>
#include <iostream>
#include <vector>
using std::cout;
using std::endl;
using std::vector;

// Constructor
Flowfield::Flowfield(Inputs inputs, MeshReader mesh_reader) {

    // Store inputs
    // TODO: Get rid of Inputs in Flowfield entirely
    this->inputs = inputs;

    // Get number of cells from the mesh
    n_cells = mesh_reader.n_cells;

    // Create initial conditions
    vector<Point> vertices;
    vector<int> vertex_ids;
    for (int i = 0; i < n_cells; i++) {

        // Vector of cell neighbor IDs
        vector<Volume*> neighbors{id_to_volume[i-1], id_to_volume[i+1]};

        // Create cell and add to map/set of cells
        // All cells are initialized to q_right for now.
        // TODO: Initialize cells with volume conditions and fix all this
        /*
        Cell *current_cell = new Cell(vertices, inputs.q_right, i, neighbors, nullptr, nullptr);
        id_to_volume[i] = current_cell;
        cells.insert(current_cell);
        volumes.insert(current_cell);
        */

    }

    // Add ghost cells
    //Ghost *left_ghost  = new Ghost(center, inputs.q_left, -1,
    //    vector<Volume*>{id_to_volume[0]}, nullptr, nullptr);
    //id_to_volume[-1] = left_ghost;
    //ghosts.insert(left_ghost);
    //volumes.insert(left_ghost);
    //center = {inputs.n_cells*inputs.dx + inputs.dx/2, .05};
    //Ghost *right_ghost = new Ghost(center, inputs.q_right, inputs.n_cells,
    //    vector<Volume*>{id_to_volume[inputs.n_cells-1]}, nullptr, nullptr);
    //id_to_volume[inputs.n_cells] = right_ghost;
    //ghosts.insert(right_ghost);
    //volumes.insert(right_ghost);

    //// Create faces of cells and assign neighbor information
    //for (int i = 0; i < inputs.n_cells + 1; i++) {
    //    // Useful points
    //    Volume *left_cell = id_to_volume[i-1];
    //    Volume *right_cell = id_to_volume[i];
    //    // Geometric information
    //    // TODO: Import this somehow instead of hardcoding
    //    vector<double> point1{i*inputs.dx, .1};
    //    vector<double> point2{i*inputs.dx, 0};
    //    Face *current_face = new Face(i, left_cell, right_cell, point1, point2);
    //    id_to_face[i] = current_face;
    //    faces.insert(current_face);
    //    // Inform neighbor cells of their new face
    //    left_cell->right_face = current_face;
    //    right_cell->left_face = current_face;
    //}

    // Start from time 0
    time = 0;

}


// Calculate fluxes through every face in domain
void Flowfield::calculate_flux(Flux &flux) {

    // For every pair of face IDs and faces, calculate the flux through each
    // of the faces
    vector<double> q_left;
    vector<double> q_right;
    vector<double> q_left_rotated;
    vector<double> q_right_rotated;
    q_left_rotated.resize(4);
    q_right_rotated.resize(4);
    for (auto &face : faces) {

        // Apply rotation matrix to the vector of conserved variables on either
        // side of the face. See Toro's book on Riemann solves, page 576, where
        // the transformation matrices are defined for 2D.
        // Rotate left state: Q^ = T*Q
        q_left = face->left_volume->q;
        q_left_rotated[0] = q_left[0];
        q_left_rotated[1] =  face->costheta * q_left[1]
            + face->sintheta * q_left[2];
        q_left_rotated[2] = -face->sintheta * q_left[1]
            + face->costheta * q_left[2];
        q_left_rotated[3] = q_left[3];
        // Rotate right state: Q^ = T*Q
        q_right = face->right_volume->q;
        q_right_rotated[0] = q_right[0];
        q_right_rotated[1] =  face->costheta * q_right[1]
            + face->sintheta * q_right[2];
        q_right_rotated[2] = -face->sintheta * q_right[1]
            + face->costheta * q_right[2];
        q_right_rotated[3] = q_right[3];

        // Get face flux using rotated left and right states. This returns the
        // rotated flux.
        vector<double> flux_rotated = flux.steger_warming(q_left_rotated, q_right_rotated,
            inputs.gamma);

        // Now, this flux has to be un-rotated. This is done by left multiplying
        // the flux vector with the inverse of the transformation matrix. Again,
        // see the Toro text mentioned above for more details of this.
        // F = (T^-1)*F^
        face->flux[0] = flux_rotated[0];
        face->flux[1] = face->costheta * flux_rotated[1]
            - face->sintheta * flux_rotated[2];
        face->flux[2] = face->sintheta * flux_rotated[1]
            + face->costheta * flux_rotated[2];
        face->flux[3] = flux_rotated[3];

    }

}


// Integrate semi-discrete equations in time
void Flowfield::apply_time_integrator() {

    // For every pair of cells and cell IDs, apply the time integrator
    for (auto &cell : cells) {

        //// Apply to every equation
        //for (int i = 0; i < inputs.n_equations; i++) {
        //    cell->q[i] = cell->q[i] - (inputs.dt/inputs.dx)
        //        * (cell->right_face->flux[i] - cell->left_face->flux[i]);
        //}
        // Update time
        time = time + inputs.dt;

    }

    // Update ghost cells
    for (auto &ghost : ghosts) {
        ghost->update();
    }

}
