#include <algorithm>
#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <math.h>
#include <vector>
#include "Flow.h"
#include "Flux.h"
using namespace Eigen;
using std::cout;
using std::endl;
using std::max_element;
using std::min_element;
using std::ofstream;
using std::vector;

Flow::Flow() {
    // Generate parameters and boundary conditions
    // TODO: Add inputfile so that inputs are not hard-coded
    n_equations = 3;
    n_cells = 200;
    n_ghosts = 2;
    cfl = .2;
    dt = 1e-4;
    time = 0;
    n_iter = 500;
    length = 1;

    // Chemistry
    gamma = 1.4;
    r_gas = 287;
}


void Flow::initialize() {
    // Allocate grid, q, f_right, and f_left
    grid.resize(n_cells + n_ghosts);
    q.resize(n_cells + n_ghosts);
    f_right.resize(n_cells + n_ghosts);
    f_left.resize(n_cells + n_ghosts);
    for (int i = 0; i < n_cells + n_ghosts; i++) {
        q[i].resize(n_equations);
        f_right[i].resize(n_equations);
        f_left[i].resize(n_equations);
    }

    // Generate grid
    for (int i=0; i<n_cells; i++) {
        grid[i] = (length*i)/(n_cells) + (length/2)/(n_cells);
    }

    // Calculate spacing
    dx = length/(n_cells);

    // Generate initial conditions
    // TODO: Add user-input initial conditions to replace hard-coded values

    // Shock
    for (int i = 0; i < n_cells/2; i++) {
        q[i](0) = 1;
        q[i](1) = 0.01;
        q[i](2) = 15000;
    }
    for (int i = n_cells/2; i < n_cells; i++) {
        q[i](0) = 1;
        q[i](1) = 0.01;
        q[i](2) = 10000;
    }

    // Calculate timestep
    //calculate_dt(q, gamma, cfl, dx);
}
/*
double Flow::calculate_dt(vector<vector<double> >& q, double gamma, double cfl, double dx) {
    // Calculate u
    vector<double> u = calculate_u(q, gamma);

    // Calculate initial time step
    double max_u = *max_element(u.begin(),u.end());
    double min_u = *min_element(u.begin(),u.end());
    double max_magnitude_u;
    if (max_u + min_u >= 0) {
        max_magnitude_u = max_u;
    } else {
        max_magnitude_u = -min_u;
    }
    double dt = cfl*dx/max_magnitude_u;
    cout << dt << endl;
    if (dt > max_dt) {
        dt = max_dt;
    }
    cout << "Timestep: " << dt << " s." << endl;
    return dt;
}
*/
void Flow::solve() {
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
        iterate(q, time, gamma, dt, dx, n_cells);
        write();
    }
}


// // // // //
// This function runs one iteration of the flow solver. This involves
// calculating wall fluxes by solving the Riemann problem, calculating the
// timestep given the CFL, and advancing the solution forward to the next
// timestep.
void Flow::iterate(
        vector<VectorXd>& q, double& time, double gamma, double cfl,
        double dx, int n_cells) {

    // Calculate timestep
    //dt = calculate_dt(q, gamma, cfl, dx);

    // Calculate fluxes, then advance every non-ghost cell forward by one
    // timestep. Indices only include non-ghost cells.
    for (int i = 1; i < n_cells + 1; i++) {
        f_right[i] = Flux::calculate_f_right(q[i], q[i+1], gamma);
        f_left[i] = Flux::calculate_f_left(q[i], q[i-1], gamma);
        q[i] = q[i] - (dt/dx)*(f_right[i] - f_left[i]);
    }

    // Update ghost cells. This is currently a reflective boundary condition.
    // TODO: Add user-input boundary condition options
    q[0] = q[1];
    q[0](1) = -q[1](1);
    q[n_cells+1] = q[n_cells];
    q[n_cells+1](1) = q[n_cells](1);

    // Advance time forward
    time = time + dt;
}


void Flow::write() {

    // Find pressure distribution
    //vector<double> p = calculate_pressure(q, gamma);

    // Find u distribution
    //vector<double> u = calculate_u(q, gamma);

    // Find temperature distribution
    //vector<double> temperature = calculate_temperature(q, gamma, r_gas);

    // Write data to solution.dat
    ofstream solution_file;
    solution_file.open("solution.dat", std::ios_base::app);
    solution_file << time << endl;
    for (int i = 0; i<n_cells; i++) {
        solution_file << grid[i] << " "
            << q[i](0) << " " << q[i](1) << " " << q[i](2) << " " << endl;
            //<< p[i] << " " << u[i] << " " << temperature[i] << endl;
    }
    solution_file << endl;
}

void Flow::output() {
    // Output results to stdout
    cout << "Final t = " << time << " s." << endl;
    cout << "Solution:" << endl;
    for (int i = 0; i < q.size(); i++) {
        cout << q[i](0) << ", " << q[i](1) << ", " << q[i](2) << "    ";
    }
    cout << endl;
}
