#include <mesh_reader.h>


// Constructor
MeshReader::MeshReader(string mesh_file) {
    this->mesh_file = mesh_file;
}


// Read mesh
void MeshReader::read_mesh() {

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
            x_coords = read_dataset<double>(file, x_coords_path);
            y_coords = read_dataset<double>(file, y_coords_path);
            connectivity = read_dataset<int>(file, connectivity_path);
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
