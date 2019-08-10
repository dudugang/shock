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
    cell_nodes.reserve(4);
    cell_faces.reserve(4);
    vector<double> q = {0, 0, 0, 0};
    for (int i = 1; i <= n_cells; i++) {

        // Find vertices of current cell from mesh connectivity
        cell_nodes = {vertices[connectivity[4*i - 4]],
                      vertices[connectivity[4*i - 3]],
                      vertices[connectivity[4*i - 2]],
                      vertices[connectivity[4*i - 1]]};

        // Create faces around boundary of current cell
        // TODO: Generalize this to work for other cell shapes
        Face *face1 = new Face(cell_nodes[0], cell_nodes[1]);
        cell_faces[0] = face1;
        Face *face2 = new Face(cell_nodes[1], cell_nodes[2]);
        cell_faces[1] = face2;
        Face *face3 = new Face(cell_nodes[2], cell_nodes[3]);
        cell_faces[2] = face3;
        Face *face4 = new Face(cell_nodes[3], cell_nodes[0]);
        cell_faces[3] = face4;

        // Add cell to mapping
        cells[i] = new Cell(cell_nodes, q, cell_faces, i);

    }


//    }

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

}


// Read mesh and BC information from HDF5 file
void MeshReader::read_hdf5() {

    // Open mesh HDF5 file
    H5File file(mesh_file, H5F_ACC_RDONLY);

    // Get base group
    string base_path = "/Base";
    Group base = file.openGroup(base_path);

    // Loop over every object inside the base group
    for (int i = 0; i < base.getNumObjs(); i++) {

        // Access object name using current iteration index
        string name = base.getObjnameByIdx(i);

        // Do nothing for \ data and Information groups
        if (name == "\ data" or name == "Information") {
            continue;
        }

        // Get mesh information from domain
        if (name == "dom-1") {
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
    dataset.read(data, PredType::NATIVE_INT, memspace, dataspace);

    // Return data
    return data;

}
