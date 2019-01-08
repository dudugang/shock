#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <math.h>
#include <vector>
#include "Flow.h"
using std::array;
using std::cout;
using std::endl;
using std::max_element;
using std::min_element;
using std::ofstream;
using std::vector;

Flow::Flow() {
    // Generate parameters and boundary conditions
    // TODO: Add inputfile so that inputs are not hard-coded
    this->n_cells = 200;
    this->cfl = .2;
    this->max_dt = .0005;
    this->time = 0;
    this->n_iter = 500;
    this->length = 1;

    // TODO: Fix this. Add ghost cells instead, this gives undefined behavior.
    // Left boundary
    this->boundary_conditions.at(0).at(0) = 1;    // rho
    this->boundary_conditions.at(1).at(0) = 0;    // rho u
    this->boundary_conditions.at(2).at(0) = 100;  // rho e
    // Right boundary
    this->boundary_conditions.at(0).at(1) = 1;    // rho
    this->boundary_conditions.at(1).at(1) = 0;    // rho u
    this->boundary_conditions.at(2).at(1) = 100;  // rho e

    // Chemistry
    this->gamma = 1.4;
    this->r_gas = 287;
}


void Flow::initialize() {
    // Allocate grid, q, q_vertex, f, and p vectors
    this->grid.resize(this->n_cells);
    this->q.resize(3);
    this->q[0].resize(this->n_cells);
    this->q[1].resize(this->n_cells);
    this->q[2].resize(this->n_cells);
    this->q_vertex.resize(3);
    this->q_vertex[0].resize(this->n_cells+1);
    this->q_vertex[1].resize(this->n_cells+1);
    this->q_vertex[2].resize(this->n_cells+1);
    this->f.resize(3);
    this->f[0].resize(this->n_cells+1);
    this->f[1].resize(this->n_cells+1);
    this->f[2].resize(this->n_cells+1);

    // Generate grid
    for (int i=0; i<this->n_cells; i++) {
        this->grid.at(i) = (this->length*i)/(n_cells) + (this->length/2)/(n_cells);
    }

    // Calculate spacing
    this->s = this->length/(this->n_cells-2);

    // Generate initial conditions
    // TODO: Add user-input initial conditions to replace hard-coded values

    /* Constant
    for (int i = 0; i<this->n_cells; i++) {
        this->q.at(0).at(i) = 1;
        this->q.at(1).at(i) = 1;
        this->q.at(2).at(i) = 100;
    } */

    // Shock
    for (int i = 0; i < this->n_cells/2; i++) {
        this->q.at(0).at(i) = 1;
        this->q.at(1).at(i) = 0;
        this->q.at(2).at(i) = 150;
    }
    for (int i = this->n_cells/2; i < this->n_cells; i++) {
        this->q.at(0).at(i) = 1;
        this->q.at(1).at(i) = 0;
        this->q.at(2).at(i) = 100;
    }

    // Square wave
    //for (int i=1; i<(this->n_cells/3); i++) {
    //    this->q.at(0).at(i) = 1.225;
    //    this->q.at(1).at(i) = 0;
    //    this->q.at(2).at(i) = 215400;
    //}
    //for (int i=this->n_cells/3; i<(2*this->n_cells/3); i++) {
    //    this->q.at(0).at(i) = 2;
    //    this->q.at(1).at(i) = 0;
    //    this->q.at(2).at(i) = 300000;
    //}
    //for (int i=2*this->n_cells/3; i<this->n_cells; i++) {
    //    this->q.at(0).at(i) = 1.225;
    //    this->q.at(1).at(i) = 0;
    //    this->q.at(2).at(i) = 215400;
    //}

    // Add boundary conditions
    this->q_vertex.at(0).at(0) = this->boundary_conditions.at(0).at(0);
    this->q_vertex.at(1).at(0) = this->boundary_conditions.at(1).at(0);
    this->q_vertex.at(2).at(0) = this->boundary_conditions.at(2).at(0);
    this->q_vertex.at(0).at(this->n_cells) = this->boundary_conditions.at(0).at(1);
    this->q_vertex.at(1).at(this->n_cells) = this->boundary_conditions.at(1).at(1);
    this->q_vertex.at(2).at(this->n_cells) = this->boundary_conditions.at(2).at(1);

    // Calculate timestep
    this->calculate_dt(this->q, this->gamma, this->cfl, this->s);
}

// // // // //
// This function takes a Q vector (rho, rho*u, rho*e) and converts it into an F
// vector, as defined by Knight (eq. 2.5)
vector<double> Flow::calculate_f_vector(vector<double>& q, double gamma) {

    // Initialize F vector
    vector<double> f;
    f.resize(3);

    // Calculate values
    f[0] = q[1];
    f[1] = q[1]*q[1]/q[0] + (gamma - 1)*(q[2] - q[1]*q[1]/(2*q[0]));
    f[2] = q[1]*q[2]/q[0] + (gamma - 1)*(q[1]/q[0])*(q[2] - q[1]*q[1]/(2*q[0]));

    return f;
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

double Flow::calculate_dt(vector<vector<double> >& q, double gamma, double cfl, double s) {
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
    double dt = cfl*s/max_magnitude_u;
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
    for (int i=1; i<=this->n_iter; i++) {
        cout << "----- Iteration " << i << ", t = " << this->time << " s. " << endl;
        this->iterate(this->q, this->q_vertex, this->f, this->time, this->gamma, this->dt, this->s, this->n_cells);
        this->write();
    }
}


// // // // //
// This function runs one iteration of the flow solver. This involves
// calculating wall fluxes by solving the Riemann problem, calculating the
// timestep given the CFL, and advancing the solution forward to the next
// timestep.
void Flow::iterate(
        vector<vector<double> >& q, vector<vector<double> >& q_vertex,
        vector<vector<double> >& f, double& time, double& gamma, double& cfl,
        double& s, int& n_cells) {

    // Calculate timestep
    dt = this->calculate_dt(this->q, this->gamma, this->cfl, this->s);

    // Calculate fluxes
    vector<vector<double> > f_right = Flux::calculate_f_right;
    vector<vector<double> > f_left = Flux::calculate_f_left;

    // Advance every grid point forward one timestep
    for (int i = 0; i < n_cells; i++) {
        q[0][i] = q[0][i] - (dt/s)*(f_right[0][i] - f_left[0][i]);
        q[1][i] = q[1][i] - (dt/s)*(f_right[1][i] - f_left[1][i]);
        q[2][i] = q[2][i] - (dt/s)*(f_right[2][i] - f_left[2][i]);
    }

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
