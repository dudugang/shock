#include <mesh_reader.h>


// Constructor
MeshReader::MeshReader(string mesh_file) {
    this->mesh_file = mesh_file;
}


// Take HDF5 mesh data and convert into data structures that the code will use
void MeshReader::create_mesh() {

    // Read data from HDF5 file
    read_hdf5();

    // Loop over and store every node in mesh
    for (int i = 1; i <= n_nodes; i++) {
        // Add to vertices mapping
        vertices[i] = Point(x_coords[i-1], y_coords[i-1]);
    }

    // Loop over every cell in mesh
    vector<Point> cell_nodes;
    vector<Face*> cell_faces;
    cell_nodes.resize(4);
    cell_faces.resize(4);
    vector<double> q = {0, 0, 0, 0}; // TODO: Initialize with ICs from hdf5 file
    for (int i = 1; i <= n_cells; i++) {

        // Find vertices of current cell from mesh connectivity
        cell_nodes = {vertices[connectivity[4*i - 4]],
                      vertices[connectivity[4*i - 3]],
                      vertices[connectivity[4*i - 2]],
                      vertices[connectivity[4*i - 1]]};

        // Create faces around boundary of current cell, and add cell to list of
        // face's neighbors
        // TODO: Generalize this to work for other cell shapes
        for (int face_index = 0; face_index < 4; face_index++) {

            // Find current face vertices. The mod operator is used so that the
            // last face loops back to the beginning of the list and uses
            // cell_nodes[3] and cell_nodes[0].
            Point vertex1 = cell_nodes[face_index % 4];
            Point vertex2 = cell_nodes[(face_index + 1) % 4];

            // Get midpoint of vertices
            Point midpoint = Geometry::find_midpoint(vertex1, vertex2);

            // Check if face already exists by looping over every face and
            // comparing midpoints
            bool face_exists = false;
            for (auto &face : faces) {
                if (face->center == midpoint) {
                    face_exists = true;
                    // If face already exists, add it to list of cell's faces
                    cell_faces[face_index] = face;
                    break;
                }
            }

            // Only create face if it doesn't already exist
            if (not face_exists) {
                // Create face
                cell_faces[face_index] = new Face(vertex1, vertex2);
            }

            // Add current cell to face's list of neighbors
            cell_faces[face_index]->neighbors.push_back(i);

        }

        // Add cell to mapping
        cells[i] = new Cell(cell_nodes, q, cell_faces, i);

        // Add faces to set
        for (auto &face : cell_faces) {
            faces.insert(face);
        }

    }

    // Create ghost cells
    int ghost_id = n_cells + 1;
    for (auto &pair : bc_connectivity) {

        // Get info from current pair in set
        string name = pair.first;
        int *conn = pair.second;

        // Loop over every face in boundary
        for (int i = 0; i < bc_face_count[name]; i++) {

            // Pointer to boundary face
            Face *boundary_face;

            // Search for existing face that has the same midpoint, then set
            // boundary face to be the same as that face, and add this ghost as
            // a neighbor to the face
            Point midpoint = Geometry::find_midpoint(vertices[conn[2*i]], vertices[conn[2*i+1]]);
            cout << "Boundary node IDs: " << conn[i] << "    " << conn[i+1] << endl;
            bool found_face = false;
            for (auto &face : faces) {
                if (face->center == midpoint) {
                    boundary_face = face;
                    boundary_face->neighbors.push_back(ghost_id);
                    found_face = true;
                    break;
                }
            }
            if (not found_face) {
                cout << "No face found to match boundary!" << endl;
                cout << "Desired midpoint: " << midpoint.x << "    " <<  midpoint.y << endl;
                cout << "Available face centers: " << endl;
                for (auto &face : faces) {
                    cout << face->center.x << "    " << face->center.y << endl;
                }
            }

            // Some random nodes for the ghost, since ghost nodes don't matter
            // TODO: Is this actually true?
            vector<Point> nodes(4, Point(0,0));
            // Put face in vector
            vector<Face*> boundary_faces = {boundary_face};
            // Create ghost cell and add to map
            ghosts[ghost_id] = new Ghost(nodes, q, boundary_faces, ghost_id);
            ghost_id++;

        }

    }

}


// Read mesh and BC information from HDF5 file
void MeshReader::read_hdf5() {

    // Open mesh HDF5 file
    H5File file(mesh_file, H5F_ACC_RDONLY);

    // Get mesh information from domain
    string x_coords_path = "/Base/dom-1/GridCoordinates/CoordinateX/\ data";
    string y_coords_path = "/Base/dom-1/GridCoordinates/CoordinateY/\ data";
    string connectivity_path = "/Base/dom-1/QuadElements/ElementConnectivity/\ data";
    string range_path = "/Base/dom-1/\ data";
    x_coords = read_dataset<double>(file, x_coords_path);
    y_coords = read_dataset<double>(file, y_coords_path);
    connectivity = read_dataset<int>(file, connectivity_path);
    int *range = read_dataset<int>(file, range_path);
    n_nodes = range[0];
    n_cells = range[1];
    // x_coords and y_coords give coordinates of nodes in mesh.
    // connectivity gives node IDs of 4 nodes per cell ID, so the first
    // 4 numbers are the 4 node IDs of cell 1, the second 4 numbers are
    // the 4 node IDs of cell 2, etc.

    // Get domain group
    string domain_path = "/Base/dom-1";
    Group domain = file.openGroup(domain_path);

    // Loop over every object inside the domain group
    for (int i = 0; i < domain.getNumObjs(); i++) {

        // Access object name using current iteration index
        string name = domain.getObjnameByIdx(i);

        // Ignore the following groups.
        unordered_set<string> ignore = {"\ data", "FamilyName", "GridCoordinates",
            "QuadElements", "ZoneBC", "ZoneType"};
        if (ignore.find(name) != ignore.end()) {
            continue;

        // Get boundary conditions
        } else {
            string bc_connectivity_path = "/Base/dom-1/" + name + "/ElementConnectivity/\ data";
            string bc_range_path = "/Base/dom-1/" + name + "/ElementRange/\ data";
            bc_connectivity[name] = read_dataset<int>(file, bc_connectivity_path);
            int *bc_range = read_dataset<int>(file, bc_range_path);
            bc_face_count[name] = bc_range[1] - bc_range[0] + 1;
        }

    }

}


// Read data from an HDF5 dataset
template<class T>
T* MeshReader::read_dataset(H5File file, string path) {

    // Get dataset and its dataspace
    DataSet dataset = file.openDataSet(path);
    DataSpace dataspace = dataset.getSpace();

    // Get rank of dataspace
    int rank = dataspace.getSimpleExtentNdims();
    // Get dimensions of dataspace
    hsize_t dimensions[rank];
    dataspace.getSimpleExtentDims(dimensions, nullptr);

    // Define the memory dataspace
    DataSpace memspace(rank, dimensions);

    // Read data into array
    T *data = new T[dimensions[0]];
    read_dataset_contents(data, dataset, memspace, dataspace);

    // Return data
    return data;

}


// Read integers from an HDF5 dataset
void MeshReader::read_dataset_contents(int *data, DataSet dataset, DataSpace memspace, DataSpace dataspace) {
    dataset.read(data, PredType::NATIVE_INT, memspace, dataspace);
}


// Read doubles from an HDF5 dataset
void MeshReader::read_dataset_contents(double *data, DataSet dataset, DataSpace memspace, DataSpace dataspace) {
    dataset.read(data, PredType::NATIVE_DOUBLE, memspace, dataspace);
}
