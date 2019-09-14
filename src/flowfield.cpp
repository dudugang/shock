#include <flowfield.h>
#include <time_integrator.h>


// Constructor
Flowfield::Flowfield(Inputs inputs, MeshReader mesh_reader) {

    // Store inputs
    // TODO: Get rid of Inputs in Flowfield entirely
    this->inputs = inputs;

    // Get number of cells/nodes from the mesh
    n_cells = mesh_reader.n_cells;
    n_nodes = mesh_reader.n_nodes;

    // Get vertices, cells, ghosts, and faces from MeshReader object
    vertices     = mesh_reader.vertices;
    connectivity = mesh_reader.connectivity;
    cells        = mesh_reader.cells;
    ghosts       = mesh_reader.ghosts;
    volumes      = mesh_reader.volumes;
    faces        = mesh_reader.faces;

    // Calculate face normal vectors
    for (auto &face : faces) {
        face->find_normal_vector(cells, ghosts, n_cells);
    }

    // Initialize ghost cell values according to boundary conditions
    for (auto &pair : ghosts) {
        // pair contains the ghost ID in pair.first and a pointer to the ghost
        // in pair.second
        pair.second->update(inputs, cells);
    }

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
    vector<double> flux_rotated;
    q_left.resize(4);
    q_right.resize(4);
    q_left_rotated.resize(4);
    q_right_rotated.resize(4);
    flux_rotated.resize(4);
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

        /*
        cout << "Face at: " << face->center.x << " " << face->center.y << endl;
        cout << "Left:" << endl;
        for (auto num : face->q_left) {
            cout << num << " ";
        }
        cout << endl;
        cout << "Right:" << endl;
        for (auto num : face->q_right) {
            cout << num << " ";
        }
        cout << endl;
        cout << "Left, rotated:" << endl;
        for (auto num : q_left_rotated) {
            cout << num << " ";
        }
        cout << endl;
        cout << "Right, rotated:" << endl;
        for (auto num : q_right_rotated) {
            cout << num << " ";
        }
        cout << endl;
        */

        // Get face flux using rotated left and right states. This returns the
        // rotated flux.
        //cout << "Face flux at: " << face->center.x << " " << face->center.y << endl;
        flux_rotated = flux.steger_warming(q_left_rotated, q_right_rotated,
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
    for (auto &face : faces) {

        // Reference to neighbor cell ID's, in increasing order
        int &smaller_id = face->neighbors[0];
        int &larger_id  = face->neighbors[1];

        // Update q on left and right of cell by setting
        face->q_left  = volumes[smaller_id]->q;
        face->q_right = volumes[larger_id]->q;

    }

}


// Integrate semi-discrete equations in time
void Flowfield::iterate(TimeIntegrator* time_integrator, Flux& flux) {

    // Integrate in time
    time_integrator->integrate(*this, cells, flux, inputs.dt);

    // Update time
    time += inputs.dt;

    // Update ghost cells to new time
    for (auto &pair : ghosts) {
        // pair contains the ghost ID in pair.first and a pointer to the ghost
        // in pair.second
        pair.second->update(inputs, cells);
    }

}
