#include <output.h>
#include <cell.h>


// Constructor
Output::Output(Inputs inputs, Flowfield &flow) {

    // Wipe old solution file
    remove("solution.dat");

    // Write header info
    ofstream solution_file;
    solution_file.open("solution.dat", std::ios_base::app);
    solution_file << "TITLE = \"Shock CFD Simulation\"" << endl;
    solution_file << "VARIABLES = \"x\", \"y\", \"Pressure\"" << endl;
    solution_file << "ZONE T=\"Iteration 0\", DATAPACKING=BLOCK, NODES="
        << flow.n_nodes << ", ELEMENTS=" << flow.n_cells
        << ", ZONETYPE=FEQUADRILATERAL, VARLOCATION=([3]=CellCentered)" << endl;

    // Write x position
    for (size_t i = 1; i <= flow.vertices.size(); i++) {
        solution_file << flow.vertices[i].x << endl;
    }
    solution_file << endl;
    // Write y position
    for (size_t i = 1; i <= flow.vertices.size(); i++) {
        solution_file << flow.vertices[i].y << endl;
    }
    solution_file << endl;
    // Write pressure
    double p;
    for (size_t i = 1; i <= flow.cells.size(); i++) {
        double rho = flow.cells[i]->q[0];
        double u = flow.cells[i]->q[1] / rho;
        double v = flow.cells[i]->q[2] / rho;
        p = (inputs.gamma-1) * (flow.cells[i]->q[3]
          - .5*rho*(u*u + v*v));
        solution_file << p << endl;
    }
    solution_file << endl;

    // Write connectivity information
    for (int cell = 0; cell < flow.n_cells; cell++) {
        // For each node in cell
        // TODO: don't hardcode number of nodes in cell
        for (int i = 0; i < 4; i++) {
            solution_file << flow.connectivity[4*cell + i] << " ";
        }
        solution_file << endl;
    }
    solution_file << endl;

}


// Add a time point to the vector of saved times to be outputted later
// TODO: I believe that means times won't be saved right if the program crashes?
void Output::add_time(double time) {

    times.push_back(time);

}


// Output information after iterations
void Output::print(Flowfield flow, int i) {

    Cell *cell = flow.cells[5];
    cout << "Iteration " << i << ": " << cell->q[0] << "  " << cell->q[1]
         << "  " << cell->q[2] << "  " << cell->q[3] << endl;

}


// Write results to file for postprocessing
void Output::write(Flowfield flow, int i) {

    // Only output every output_rate iterations
    if (i % flow.inputs.output_rate == 0) {

        // Write data to solution.dat
        ofstream solution_file;
        solution_file.open("solution.dat", std::ios_base::app);

        // Write header information
        solution_file << "ZONE T=\"Iteration " << i+1 << "\", DATAPACKING=BLOCK, NODES="
            << flow.n_nodes << ", ELEMENTS=" << flow.n_cells
            << ", ZONETYPE=FEQUADRILATERAL, VARLOCATION=([3]=CellCentered), "
            << "VARSHARELIST=([1-2]=1), CONNECTIVITYSHAREZONE=1, "
            << "SOLUTIONTIME=" << flow.time << endl;

        // Write pressure
        double p;
        for (size_t i = 1; i <= flow.cells.size(); i++) {
            double rho = flow.cells[i]->q[0];
            double u = flow.cells[i]->q[1] / rho;
            double v = flow.cells[i]->q[2] / rho;
            p = (flow.inputs.gamma-1) * (flow.cells[i]->q[3]
              - .5*rho*(u*u + v*v));
            solution_file << p << endl;
        }
        solution_file << endl << endl;

    }

}


// Write final results to HDF5 dataset
void Output::write_results(string case_file, unordered_map<int, Cell*> &cells, unsigned int n_cells) {

    // Get data
    double rho[n_cells];
    for (auto &pair : cells) {

        // Convenient variables
        int id = pair.first;
        Cell *cell = pair.second;

        // Extract results
        rho[id - 1] = cell->q[0];

    }

    // Open HDF5 file
    H5File file(case_file, H5F_ACC_RDWR);
    // Create groups
    //Group base_iterative_data(file.createGroup("/Base/BaseIterativeData"));
    //Group zone1(file.createGroup("/Base/Zone1"));
    Group flow_solution(file.createGroup("/Base/dom-1/FlowSolution"));
    Group density(file.createGroup("/Base/dom-1/FlowSolution/Density"));

    // Write all data to file
    string times_path = "/Base/BaseIterativeData/TimeValues";
    string n_steps_path = "/Base/BaseIterativeData/NumberOfSteps";
    string rho_path = "/Base/dom-1/FlowSolution/Density/ data";
    //unsigned int n_steps[] = {times.size()};
    //write_dataset(file, times_path, &times[0], times.size());
    //write_dataset(file, n_steps_path, n_steps, 1);
    //write_dataset(file, rho_path, rho, n_cells);
    //file.link(H5G_LINK_HARD, "Base/dom-1/FamilyName", "Base/Zone1/FamilyName");
    //file.link(H5G_LINK_HARD, "Base/dom-1/ZoneBC", "Base/Zone1/ZoneBC");
    //file.link(H5G_LINK_HARD, "Base/dom-1/ZoneType", "Base/Zone1/ZoneType");
    //file.link(H5G_LINK_HARD, "Base/dom-1/GridCoordinates", "Base/Zone1/GridCoordinates");
    //file.link(H5G_LINK_HARD, "Base/dom-1/QuadElements", "Base/Zone1/QuadElements");
    //file.link(H5G_LINK_HARD, "Base/dom-1/ data", "Base/Zone1/ data");
    //file.link(H5G_LINK_HARD, "Base/dom-1/in", "Base/Zone1/in");
    //file.link(H5G_LINK_HARD, "Base/dom-1/out", "Base/Zone1/out");
    //file.link(H5G_LINK_HARD, "Base/dom-1/sym", "Base/Zone1/sym");
    //file.link(H5G_LINK_HARD, "Base/dom-1/wedge", "Base/Zone1/wedge");

}


// Write integer data to an HDF5 dataset
void Output::write_dataset(H5File file, string dataset_path, int data[], hsize_t length) {

    // Create dataspace
    int rank = 1;
    hsize_t dims[] = {length};
    DataSpace dataspace(rank, dims);

    // Create dataset
    // TODO: Might need an integer and double version of this
    DataSet dataset = file.createDataSet(dataset_path, PredType::NATIVE_INT, dataspace);

    // Write data to dataset
    // TODO: Might need an integer and double version of this, too
    dataset.write(data, PredType::NATIVE_INT, dataspace);

}


// Write double precision data to an HDF5 dataset
void Output::write_dataset(H5File file, string dataset_path, double data[], hsize_t length) {

    // Create dataspace
    int rank = 1;
    hsize_t dims[] = {length};
    DataSpace dataspace(rank, dims);

    // Create dataset
    // TODO: Might need an integer and double version of this
    DataSet dataset = file.createDataSet(dataset_path, PredType::NATIVE_DOUBLE, dataspace);

    // Write data to dataset
    // TODO: Might need an integer and double version of this, too
    dataset.write(data, PredType::NATIVE_DOUBLE, dataspace);

}
