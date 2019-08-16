#include <output.h>
#include <cell.h>


// Constructor
Output::Output(Inputs inputs, int n_cells) {

    // Wipe old solution file
    remove("solution.dat");

    // Calculate total number of file writes: n_iterations / output_rate,
    // rounded up, plus one (since initial conditions always written).
    // Note: the following is basically just a hacky way of doing ceiling
    // integer division (i.e. rounding up).
    int n_writes = (inputs.n_iterations + inputs.output_rate - 1) / inputs.output_rate + 1;

    // Write header info
    ofstream solution_file;
    solution_file.open("solution.dat", std::ios_base::app);
    solution_file << n_cells << " " << n_writes << endl;

}


// Add a time point to the vector of saved times to be outputted later
// TODO: I believe that means times won't be saved right if the program crashes?
void Output::add_time(double time) {

    times.push_back(time);

}


// Output information after iterations
void Output::print(Flowfield flow, int i) {

    //cout << "Iteration " << i << ": " << flow.cells[5]->q[0] << endl;

}


// Output results to stdout
void Output::final_print(Flowfield flow) {

    cout << "Final t = " << flow.time << " s." << endl;
    cout << "Solution:" << endl;
    for (auto &pair : flow.cells) {
        Cell* cell = pair.second;
        cout << "Volume ID " << cell->id << ": " << cell->q[0] << ", "
             << cell->q[1] << ", " << cell->q[2] << endl;
    }
    cout << endl;

}

// Write results to file for postprocessing
void Output::write(Flowfield flow, int i) {

    // Only output every output_rate iterations
    if (i % flow.inputs.output_rate == 0) {
        // Write data to solution.dat
        ofstream solution_file;
        solution_file.open("solution.dat", std::ios_base::app);
        solution_file << flow.time << endl;
        for (auto &pair : flow.cells) {
            Cell* cell = pair.second;
            solution_file << cell->center.x << " " << cell->center.y << " "
                << cell->q[0] << " " << cell->q[1] << " " << cell->q[2] << " "
                << cell->q[3] << endl;
        }
        solution_file << endl;
    }

}


// Write final results to HDF5 dataset
void Output::write_results(string case_file, unordered_map<int, Cell*> &cells, int n_cells) {

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
    Group base_iterative_data(file.createGroup("/Base/BaseIterativeData"));

    // Write all data to file
    string rho_path = "/Base/BaseIterativeData";
    string times_path = "/Base/BaseIterativeData/TimeValues";
    string n_steps_path = "/Base/BaseIterativeData/NumberOfSteps";
    int n_steps[] = {times.size()};
    //write_dataset(file, rho_path, rho, n_cells);
    write_dataset(file, times_path, &times[0], times.size());
    write_dataset(file, n_steps_path, n_steps, 1);

}


// Write integer data to an HDF5 dataset
void Output::write_dataset(H5File file, string dataset_path, int data[], int length) {

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
void Output::write_dataset(H5File file, string dataset_path, double data[], int length) {

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
