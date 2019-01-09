#include <algorithm>
#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <math.h>
#include <vector>
#include "Flow.h"
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
    this->n_equations = 3;
    this->n_cells = 200;
    this->n_ghosts = 2;
    this->cfl = .2;
    this->max_dt = .0005;
    this->time = 0;
    this->n_iter = 500;
    this->length = 1;

    // Chemistry
    this->gamma = 1.4;
    this->r_gas = 287;
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
        grid.at(i) = (length*i)/(n_cells) + (length/2)/(n_cells);
    }

    // Calculate spacing
    dx = length/(n_cells);

    // Generate initial conditions
    // TODO: Add user-input initial conditions to replace hard-coded values

    // Shock
    for (int i = 0; i < n_cells/2; i++) {
        q.at(0).at(i) = 1;
        q.at(1).at(i) = 0;
        q.at(2).at(i) = 150;
    }
    for (int i = n_cells/2; i < n_cells; i++) {
        q.at(0).at(i) = 1;
        q.at(1).at(i) = 0;
        q.at(2).at(i) = 100;
    }

    // Calculate timestep
    calculate_dt(q, gamma, cfl, dx);
}

vector<double> Flow::calculate_pressure(vector<vector<double> >& q, double gamma) {
    // Calculate pressure from a given Q vector

    // Initialize pressure vector
    vector<double> p;
    p.resize(q[0].size());

    // Calculate values
    for (int i=0; i<q[0].size(); i++) {
        p[i] = (gamma - 1)*(q.at(2).at(i) - q.at(1).at(i)*q.at(1).at(i)/(2*q.at(0).at(i)));
    }

    return p;
}

vector<double> Flow::calculate_u(vector<vector<double> >& q, double gamma) {
    // Calculate u (x-velocity) from a given Q vector

    // Initialize u vector
    vector<double> u;
    u.resize(q[0].size());

    // Calculate values
    for (int i=0; i<q[0].size(); i++) {
        u[i] = q.at(1).at(i) / q.at(0).at(i);
    }

    return u;
}

vector<double> Flow::calculate_temperature(vector<vector<double> >& q, double gamma, double r_gas) {
    // Calculate temperature from a given Q vector

    // Initialize temperature vector
    vector<double> temperature;
    temperature.resize(q[0].size());

    // Calculate values
    for (int i=0; i<q[0].size(); i++) {
        temperature[i] = ((gamma-1)/r_gas) * (q.at(2).at(i)/q.at(0).at(i));
    }

    return temperature;
}

double Flow::calculate_dt(vector<vector<double> >& q, double gamma, double cfl, double dx) {
    // Calculate u
    vector<double> u = this->calculate_u(q, gamma);

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

void Flow::solve() {
    // Wipe old solution file
    remove("solution.dat");

    // Write header info
    ofstream solution_file;
    solution_file.open("solution.dat", std::ios_base::app);
    solution_file << this->n_cells << " " << this->n_iter << endl;

    // Write initial conditions
    this->write();

    // Iterate until n_iter is reached
    for (int i=1; i<=n_iter; i++) {
        cout << "----- Iteration " << i << ", t = " << time << " s. " << endl;
        iterate(q, q_vertex, f, time, gamma, dt, dx, n_cells);
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
    dt = this->calculate_dt(q, gamma, cfl, dx);

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
    vector<double> p = this->calculate_pressure(this->q, this->gamma);

    // Find u distribution
    vector<double> u = this->calculate_u(this->q, this->gamma);

    // Find temperature distribution
    vector<double> temperature = this->calculate_temperature(this->q, this->gamma, this->r_gas);

    // Write data to solution.dat
    ofstream solution_file;
    solution_file.open("solution.dat", std::ios_base::app);
    solution_file << this->time << endl;
    for (int i = 0; i<this->n_cells; i++) {
        solution_file << this->grid[i] << " "
            << this->q.at(0).at(i) << " " << this->q.at(1).at(i) << " " << this->q.at(2).at(i) << " "
            << p[i] << " " << u[i] << " " << temperature[i] << endl;
    }
    solution_file << endl;
}

void Flow::output() {
    // Output results to stdout
    cout << "Final t = " << this->time << " s." << endl;
    cout << "Solution:" << endl;
    for (int i = 0; i < this->q[0].size(); i++) {
        cout << this->q.at(0).at(i) << ", " << this->q.at(1).at(i) << ", " << this->q.at(2).at(i) << "    ";
    }
    cout << endl;
}
