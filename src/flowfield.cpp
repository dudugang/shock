#include <Flow.h>
#include <Flux.h>
#include <algorithm>
#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <math.h>
#include <vector>
using namespace Eigen;
using std::cout;
using std::endl;
using std::max_element;
using std::min_element;
using std::ofstream;
using std::vector;

Flowfield::Flowfield() {
    // Generate parameters and boundary conditions
    // TODO: Add inputfile so that inputs are not hard-coded
    n_equations = 3;
    n_cells = 200;
    n_ghosts = 2;
    cfl = .2;
    dt = 1e-6;
    time = 0;
    n_iter = 2000;
    length = 1;

    // Chemistry
    gamma = 1.4;
    r_gas = 287;
}


void Flowfield::initialize() {
    // Allocate space for data
    grid.resize(n_cells + n_ghosts);
    q.resize(n_cells + n_ghosts);
    q_old.resize(n_cells + n_ghosts);
    f_right.resize(n_cells + n_ghosts);
    f_left.resize(n_cells + n_ghosts);
    for (int i = 0; i < n_cells + n_ghosts; i++) {
        q[i].resize(n_equations);
        q_old[i].resize(n_equations);
        f_right[i].resize(n_equations);
        f_left[i].resize(n_equations);
    }

    // Generate grid
    for (int i=0; i < (n_cells + n_ghosts); i++) {
        grid[i] = (length*i)/(n_cells) - (length/2)/(n_cells);
    }

    // Calculate spacing
    dx = length/(n_cells);

    // Generate initial conditions
    // TODO: Add user-input initial conditions to replace hard-coded values

    // Shock
    for (int i = 0; i < (n_cells + n_ghosts)/2; i++) {
        q[i](0) = 1.5;
        q[i](1) = 0.00001;
        q[i](2) = 3.7e6;
    }
    for (int i = (n_cells + n_ghosts)/2; i < (n_cells + n_ghosts); i++) {
        q[i](0) = 1;
        q[i](1) = 0.00001;
        q[i](2) = 2.5e6;
    }

    // Calculate timestep
    //calculate_dt(q, gamma, cfl, dx);
}


void Flowfield::solve() {
    // Wipe old solution file
    remove("solution.dat");

    // Write header info
    ofstream solution_file;
    solution_file.open("solution.dat", std::ios_base::app);
    solution_file << n_cells << " " << n_iter << endl;

    // Write initial conditions
    write();

    // Iterate until n_iter is reached
    for (int i=1; i<=n_iter; i++) {
        cout << "----- Iteration " << i << ", t = " << time << " s. " << endl;
        iterate();
        write();
    }
}


void Flowfield::write() {

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

void Flowfield::output() {
    // Output results to stdout
    cout << "Final t = " << time << " s." << endl;
    cout << "Solution:" << endl;
    for (int i = 0; i < (n_cells + n_ghosts); i++) {
        cout << q[i](0) << ", " << q[i](1) << ", " << q[i](2) << "    ";
    }
    cout << endl;
}
