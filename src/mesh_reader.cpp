#include <mesh_reader.h>
#include <algorithm>
#include <cell.h>
#include <face.h>
#include <ghost.h>
#include <state.h>
#include <volume.h>


// Constructor
MeshReader::MeshReader(string mesh_file) {
    this->mesh_file = mesh_file;
}


// Take HDF5 mesh data and convert into data structures that the code will use
void MeshReader::create_mesh(Inputs &inputs) {

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

    // TODO: Initialize with ICs from hdf5 file
    vector<double> q;
    q.resize(4);
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

            // Create face and add to list of cell's faces
            cell_faces[face_index] = new Face(vertex1, vertex2);

            // Add current cell to face's list of neighbors
            cell_faces[face_index]->neighbors.push_back(i);

        }

        // Add cell to mapping
        // TODO: Fix this hack
        if (cell_nodes[0].x < .5) {
            q = State::physical_to_conserved(inputs.vc["driver"], inputs.gamma);
        } else {
            q = State::physical_to_conserved(inputs.vc["driven"], inputs.gamma);
        }
        cells[i] = new Cell(cell_nodes, q, cell_faces, i);
        volumes[i] = cells[i];

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
            for (auto &face : faces) {
                if (face->center == midpoint) {
                    boundary_face = face;
                    boundary_face->neighbors.push_back(ghost_id);
                    break;
                }
            }

            // Ghost only has two nodes, which are the nodes if its boundary
            // face
            vector<Point> nodes = {boundary_face->point1, boundary_face->point2};
            // Put face in vector
            vector<Face*> boundary_faces = {boundary_face};
            // Create ghost cell, add to map, and save boundary condition name
            ghosts[ghost_id] = new Ghost(nodes, q, boundary_faces, ghost_id);
            volumes[ghost_id] = ghosts[ghost_id];
            ghosts[ghost_id]->type = name;

            ghost_id++;

        }

    }

    // Find neighbors of every cell
    find_cell_neighbors();

    // Find neighbors of every face
    find_face_neighbors();

    // Sort each face's neighbors in ascending order of volume ID's
    for (auto &face : faces) {
        face->sort_neighbors();
    }

    // Combine duplicate faces using connectivity information
    combine_duplicate_faces();

}


// Finds the neighbors of every cell/ghost in the mesh
void MeshReader::find_cell_neighbors() {

    // Vector of all volumes
    vector<Volume*> volumes;
    volumes.reserve(cells.size() + ghosts.size());
    // Put all cells in vector of volumes
    for (auto &pair : cells)  { volumes.push_back(pair.second); }
    // Put all ghosts in vector of volumes
    for (auto &pair : ghosts) { volumes.push_back(pair.second); }

    // Sort all volumes by the x-coordinate of their first face's center
    std::sort(volumes.begin(), volumes.end(), [](Volume* v1, Volume* v2) {
        return v1->faces[0]->center.x < v2->faces[0]->center.x;
    });

    // Loop through all volumes
    for (int current_index = 0; current_index < volumes.size(); current_index++) {

        // Convenient point to current volume
        Volume *current = volumes[current_index];

        // Loop through volumes again, looking for the neighbors of each cell.
        // However, since cells are sorted by one of their face's centers,
        // looking for cells near the same index is most likely to find the
        // match. This should save lots of time.
        for (int j = 1; j < volumes.size(); j++) {

            // Convert j to id of candidate cell, to be determined if this is a
            // neighbor
            int candidate_index;
            if (j % 2 != 0) {
                candidate_index = current_index + (j + 1)/2;
            } else {
                candidate_index = current_index - j/2;
            }

            // If candidate index is negative, skip this iteration
            if (candidate_index < 0) { continue; }
            // If candidate index is the same as the current cell index, skip
            // this iteration
            if (candidate_index == current_index) { continue; }
            // If candidate index is greater than or equal to the total number
            // of volumes, skip this iteration
            if (candidate_index >= volumes.size()) { continue; }

            // Convenient pointer to candidate volume
            Volume *candidate = volumes[candidate_index];

            // Loop over every vertex in the candidate cell
            int matches = 0;
            for (auto &vertex : candidate->vertices) {

                // If this vertex of the candidate cell matches a vertex in the
                // current cell, then increase matches
                if (std::find(current->vertices.begin(),
                    current->vertices.end(), vertex) != current->vertices.end()) {
                    matches++;
                }

            }

            // If there are two matches, then there are two vertices shared
            // between the two volumes, which means that they are neighbors.
            // This is only true in 2D.
            if (matches == 2) { current->neighbors.push_back(candidate->id); }

            // Stop looking for a volume's neighbors if all neighbors have already
            // been found
            // TODO: Make this work for any shape of cell
            if (current->type == "flow"
                and current->neighbors.size() == 4) {
                break;
            } else if (current->type == "ghost"
                and current->neighbors.size() == 1) {
                break;
            }

        }

    }

}


// Find the other cell neighboring each face.
void MeshReader::find_face_neighbors() {

    // Loop over every face
    for (auto &face : faces) {

        // Pointer to cell adjacent to face
        Cell *cell = cells[face->neighbors[0]];

        // Loop over every neighboring cell
        bool found_neighbor = false;
        for (auto &neighbor_id : cell->neighbors) {

            // Pointer to the neighbor volume
            Volume *neighbor;
            if (neighbor_id <= n_cells) {
                neighbor = cells[neighbor_id];
            } else {
                neighbor = ghosts[neighbor_id];
            }

            // Loop over every face of the neighbor
            for (auto &neighbor_face : neighbor->faces) {
                // If the neighbor's face matches the original face, then add
                // the neighbor to the original face's list of neighbors
                if (*face == *neighbor_face) {
                    face->neighbors.push_back(neighbor->id);
                    found_neighbor = true;
                    break;
                }
            }

            // Stop this if the neighbor was found
            if (found_neighbor) { break; }

        }

    }

}


// Find all face duplicates created during initialization of volumes and get rid
// of the extra, replacing it with a pointer to the original. This way, the flux
// calculations don't happen twice.
void MeshReader::combine_duplicate_faces() {

    // Loop over every volume
    for (auto &pair : volumes) {

        // Cell id and pointer to cell
        int volume_id = pair.first;
        Volume* volume = pair.second;

        // Loop over all faces of volume
        for (int i = 0; i < volume->faces.size(); i++) {

            // Pointer to this face
            Face* face = volume->faces[i];

            // If this cell is the second neighbor of the face, then it has
            // higher cell ID than the adjacent cell, so do nothing for this
            // face
            int adjacent_id;
            if (face->neighbors[1] == volume_id or volume_id > n_cells) {
                continue;
            // Otherwise, this cell is the first neighbor, meaning it has lower
            // cell ID than its adjacent cell.
            } else {
                adjacent_id = face->neighbors[1];
            }

            // Get pointer to adjacent cell
            Volume *adjacent = volumes[adjacent_id];

            // If adjacent is a ghost, do nothing
            if (adjacent->type != "flow") { continue; }

            // Loop over faces of adjacent cell
            for (int j = 0; j < adjacent->faces.size(); j++) {
                // If face from original cell matches this adjacent face
                if (*face == *(adjacent->faces[j])) {
                    // Delete original face
                    delete face;
                    // Set pointer equal to new face
                    volume->faces[i] = adjacent->faces[j];
                }
            }

        }

    }

    // Clear faces from set
    faces.clear();
    // Loop over every cell
    for (auto &pair : cells) {
        // Cell pointer and its ID
        Cell *cell = pair.second;
        int id = pair.first;
        // Loop over every face
        for (auto &face : cell->faces) {
            // If this cell's ID is less than the adjacent cell, then add the
            // face to the set
            if (face->neighbors[0] == id) {
                faces.insert(face);
            }
        }
    }


}


// Read mesh and BC information from HDF5 file
void MeshReader::read_hdf5() {

    // Open mesh HDF5 file
    H5File file(mesh_file, H5F_ACC_RDONLY);

    // Get mesh information from domain
    string x_coords_path = "/Base/dom-1/GridCoordinates/CoordinateX/ data";
    string y_coords_path = "/Base/dom-1/GridCoordinates/CoordinateY/ data";
    string range_path = "/Base/dom-1/ data";
    string connectivity_path = "/Base/dom-1/QuadElements/ElementConnectivity/ data";
    x_coords = read_dataset<double>(file, x_coords_path);
    y_coords = read_dataset<double>(file, y_coords_path);
    int *range = read_dataset<int>(file, range_path);
    n_nodes = range[0];
    n_cells = range[1];
    conn = read_dataset<int>(file, connectivity_path);
    connectivity.resize(n_cells*4);
    for (int i = 0; i < n_cells*4; i++) {
        connectivity[i] = conn[i];
    }
    // x_coords and y_coords give coordinates of nodes in mesh.
    // connectivity gives node IDs of 4 nodes per cell ID, so the first
    // 4 numbers are the 4 node IDs of cell 1, the second 4 numbers are
    // the 4 node IDs of cell 2, etc.

    // Get domain group
    string domain_path = "/Base/dom-1";
    Group domain = file.openGroup(domain_path);

    // Loop over every object inside the domain group
    for (unsigned int i = 0; i < domain.getNumObjs(); i++) {

        // Access object name using current iteration index
        string name = domain.getObjnameByIdx(i);

        // Ignore the following groups.
        unordered_set<string> ignore = {" data", "FamilyName", "GridCoordinates",
            "QuadElements", "ZoneBC", "ZoneType"};
        if (ignore.find(name) != ignore.end()) {
            continue;

        // Get boundary conditions
        } else {
            string bc_connectivity_path = "/Base/dom-1/" + name + "/ElementConnectivity/ data";
            string bc_range_path = "/Base/dom-1/" + name + "/ElementRange/ data";
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
