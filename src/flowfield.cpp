#include <flowfield.h>


// Constructor
Flowfield::Flowfield(Inputs inputs, MeshReader mesh_reader) {

    // Store inputs
    // TODO: Get rid of Inputs in Flowfield entirely
    this->inputs = inputs;

    // Get number of cells/nodes from the mesh
    n_cells = mesh_reader.n_cells;

    // Get cells, ghosts, and faces from MeshReader object
    cells  = mesh_reader.cells;
    ghosts = mesh_reader.ghosts;
    faces  = mesh_reader.faces;

    // Start from time 0
    time = 0;

}


// Calculate fluxes through every face in domain
void Flowfield::calculate_flux(Flux &flux) {

    // Initialize vectors for holding left and right states of each face
    vector<double> q_left;
    vector<double> q_right;
    vector<double> q_left_rotated;
    vector<double> q_right_rotated;
    q_left_rotated.resize(4);
    q_right_rotated.resize(4);
    // Loop over all faces
    for (auto &face : faces) {

        // Apply rotation matrix to the vector of conserved variables on either
        // side of the face. See Toro's book on Riemann solves, page 576, where
        // the transformation matrices are defined for 2D.
        // Rotate left state: Q^ = T*Q
        q_left = face->q_left;
        q_left_rotated[0] = q_left[0];
        q_left_rotated[1] =  face->costheta * q_left[1]
            + face->sintheta * q_left[2];
        q_left_rotated[2] = -face->sintheta * q_left[1]
            + face->costheta * q_left[2];
        q_left_rotated[3] = q_left[3];

        // Rotate right state: Q^ = T*Q
        q_right = face->q_right;
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


// Reconstruct flow from cell centers to faces
void Flowfield::apply_reconstruction() {

    // Assume 1D reconstruction
    // TODO: Stop assuming 1D reconstruction
    // Loop over every cell
    for (auto &pair : cells) {

        // Define convenient pointer
        Cell* cell = pair.second;

        // Loop over every face in the cell
        for (auto &face : cell->faces) {

            // Find angle from horizontal between cell center and face point
            double angle = std::atan2(face->center.y - cell->center.y,
                face->center.x - cell->center.x);

            // Check if face normal vector points into cell or out of cell
            if (std::abs(angle - face->theta) < 90) {

            }

        }

    }

}


// Integrate semi-discrete equations in time
void Flowfield::apply_time_integrator() {

    // Loop over every pair of cell IDs and cells
    for (auto &pair : cells) {

        // Define convenient pointer
        Cell* cell = pair.second;

        // Apply to every equation
        double flux_integral;
        double coefficient;
        for (int i = 0; i < inputs.n_equations; i++) {
            coefficient = - (inputs.dt/cell->volume);
            // Integrate flux over every cell face
            for (auto &face : cell->faces) {
                flux_integral += face->flux[i] * face->area;
            }
            cell->q[i] += coefficient * flux_integral;
        }

        // Update time
        time = time + inputs.dt;

    }

    // Update ghost cells
    for (auto &ghost : ghosts) {
        //ghost->update();
    }

}
