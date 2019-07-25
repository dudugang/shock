#include <output.h>
using std::cout;
using std::endl;

// Output results to stdout
void Output::print(Flowfield flow) {

    cout << "Final t = " << flow.time << " s." << endl;
    cout << "Solution:" << endl;
    for (auto &pair : flow.cells) {
        Cell &cell = *(pair.second);
        cout << "Cell ID " << cell.cell_id << ": " << cell.q[0] << ", "
             << cell.q[1] << ", " << cell.q[2] << endl;
    }
    cout << endl;

}

/*
void Output::write() {

    // Find pressure distribution
    vector<double> p;
    p.resize(n_cells + n_ghosts);
    for (int i = 0; i < (n_cells + n_ghosts); i++) {
        p[i] = calculate_pressure(q[i], gamma);
    }

    // Find u distribution
    //vector<double> u = calculate_u(q, gamma);

    // Find temperature distribution
    //vector<double> temperature = calculate_temperature(q, gamma, r_gas);

    // Write data to solution.dat
    ofstream solution_file;
    solution_file.open("solution.dat", std::ios_base::app);
    solution_file << time << endl;
    for (int i = 0; i < (n_cells + n_ghosts); i++) {
        solution_file << grid[i] << " "
            << q[i](0) << " " << q[i](1) << " " << q[i](2) << " " << p[i] << endl;
            //u[i] << " " << temperature[i] << endl;
    }
    solution_file << endl;
}
*/
