#include <inputs.h>

// Constructor. This stores user inputs to the code.
Inputs::Inputs() {

    // Parameters
    dt = .000001;
    n_cells = 200;
    n_iterations = 2000;
    length = 1;
    dx = length/n_cells;
    gamma = 1.4;
    n_equations = 4;
    output_rate = 10;

    // Initial conditions of conserved variables rho, rho*u, and rho*e
    q_left.resize(4);
    q_right.resize(4);
    q_left[0]  = 1.5;
    q_left[1]  = 0.00001;
    q_left[2]  = 0.00001;
    q_left[3]  = 3.7e6;
    q_right[0] = 1;
    q_right[1] = 0.00001;
    q_right[2] = 0.00001;
    q_right[3] = 2.5e6;
}


// Read in mesh
void Inputs::read_mesh(){

    string file_name = "wedge_30deg.cgns";
    string dset = "\ format";
    int buf[15];
    hsize_t dimsm[1] = {15};
    int RANK_OUT = 1;

    H5File file(file_name, H5F_ACC_RDONLY);
    DataSet dataset = file.openDataSet(dset);
    DataSpace dataspace = dataset.getSpace();
    DataSpace memspace( RANK_OUT, dimsm );
    dataset.read(buf, PredType::NATIVE_INT, memspace, dataspace);
    cout << buf[0] << endl;

}
