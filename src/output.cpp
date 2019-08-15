#include <output.h>


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
