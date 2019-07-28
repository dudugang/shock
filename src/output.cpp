#include <output.h>
using std::cout;
using std::endl;
using std::ofstream;


// Constructor
Output::Output(Inputs inputs) {

    // Wipe old solution file
    remove("solution.dat");

    // Write header info
    ofstream solution_file;
    solution_file.open("solution.dat", std::ios_base::app);
    solution_file << inputs.n_cells << " " << inputs.n_iterations << endl;

}


// Output information after each iteration
void Output::print(Flowfield flow, int i) {

    cout << "Iteration " << i << ": " << flow.id_to_volume[5]->q[2] << endl;

}


// Output results to stdout
void Output::final_print(Flowfield flow) {

    cout << "Final t = " << flow.time << " s." << endl;
    cout << "Solution:" << endl;
    for (auto &cell : flow.cells) {
        cout << "Volume ID " << cell->volume_id << ": " << cell->q[0] << ", "
             << cell->q[1] << ", " << cell->q[2] << endl;
    }
    cout << endl;

}

// Write results to file for postprocessing
void Output::write(Flowfield flow) {

    // Find pressure distribution
    //vector<double> p;
    //p.resize(n_cells + n_ghosts);
    //for (int i = 0; i < (n_cells + n_ghosts); i++) {
    //    p[i] = calculate_pressure(q[i], gamma);
    //}
    // Find u distribution
    //vector<double> u = calculate_u(q, gamma);
    // Find temperature distribution
    //vector<double> temperature = calculate_temperature(q, gamma, r_gas);

    // Write data to solution.dat
    ofstream solution_file;
    solution_file.open("solution.dat", std::ios_base::app);
    solution_file << flow.time << endl;
    for (auto &volume : flow.volumes) {
        solution_file << volume->x << " "
            << volume->q[0] << " " << volume->q[1] << " " << volume->q[2] << " "
            << endl;
            //u[i] << " " << temperature[i] << endl;
    }
    solution_file << endl;
}
