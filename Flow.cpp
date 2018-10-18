#include <array>
#include <iostream>
#include <fstream>
#include <vector>
#include "Flow.h"
using std::array;
using std::cout;
using std::endl;
using std::ofstream;
using std::vector;

Flow::Flow() {

    // Generate parameters and boundary conditions
    // TODO: Add inputfile so that inputs are not hard-coded
    this->n_cells = 4;
    this->dt = .005;
    this->time = 0;
    this->n_iter = 1000;
    this->length = 10;
    // Left boundary
    this->boundary_conditions.at(0).at(0) = 1.225;    // rho
    this->boundary_conditions.at(1).at(0) = 0;       // rho u
    this->boundary_conditions.at(2).at(0) = 300.000;   // rho e
    // Right boundary
    this->boundary_conditions.at(0).at(1) = 1.225;    // rho
    this->boundary_conditions.at(1).at(1) = 0;       // rho u
    this->boundary_conditions.at(2).at(1) = 215.400;   // rho e
    this->gamma = 1.4;

}


void Flow::initialize() {

    // Allocate grid, q, f, and p vectors
    this->grid.resize(this->n_cells);
    this->q.resize(3);
    this->q[0].resize(this->n_cells);
    this->q[1].resize(this->n_cells);
    this->q[2].resize(this->n_cells);
    this->f.resize(3);
    this->f[0].resize(this->n_cells);
    this->f[1].resize(this->n_cells);
    this->f[2].resize(this->n_cells);
    this->p.resize(this->n_cells);

    // Generate grid
    for (int i=0; i<this->n_cells; i++) {
        this->grid.at(i) = (this->length*i)/(n_cells-2) - (this->length/2)/(n_cells-2);
    }
    for (int i = 0; i < this->grid.size(); i++) {
        cout << grid[i] << "    ";
    }
    // Calculate spacing
    this->s = this->length/(this->n_cells-2);

    // Generate initial conditions
    // TODO: Add user-input initial conditions to replace hard-coded values

    // Constant
    for (int i = 1; i<(this->n_cells)-1; i++) {
        this->q.at(0).at(i) = 1.225;
        this->q.at(1).at(i) = 0;
        this->q.at(2).at(i) = 215.400;
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
    this->q.at(0).at(0) = this->boundary_conditions.at(0).at(0);
    this->q.at(1).at(0) = this->boundary_conditions.at(1).at(0);
    this->q.at(2).at(0) = this->boundary_conditions.at(2).at(0);
    this->q.at(0).at(this->n_cells-1) = this->boundary_conditions.at(0).at(1);
    this->q.at(1).at(this->n_cells-1) = this->boundary_conditions.at(1).at(1);
    this->q.at(2).at(this->n_cells-1) = this->boundary_conditions.at(2).at(1);

}

vector<vector<double> > Flow::calculate_f_vector(vector<vector<double> >& q, int n_cells, double gamma) {
    // Calculate F vector from a given Q vector

    // Initialize F vector
    vector<vector<double> > f;
    f.resize(3);
    f[0].resize(n_cells);
    f[1].resize(n_cells);
    f[2].resize(n_cells);

    // Calculate values
    for (int i=0; i<q[0].size(); i++) {
        f.at(0).at(i) = q.at(1).at(i);
        f.at(1).at(i) = q.at(1).at(i)*q.at(1).at(i)/q.at(0).at(i) + (gamma - 1)*(q.at(2).at(i) - q.at(1).at(i)*q.at(1).at(i)/(2*q.at(0).at(i)));
        f.at(2).at(i) = q.at(1).at(i)*q.at(2).at(i)/q.at(0).at(i) + (gamma - 1)*(q.at(1).at(i)/q.at(0).at(i))*(q.at(2).at(i) - q.at(1).at(i)*q.at(1).at(i)/(2*q.at(0).at(i)));
    }
    cout << f[0][1] << "    " << f[1][1] << "    " << f[2][1] << endl;
    this->f = f;
    return f;

}

vector<double> Flow::calculate_pressure(vector<vector<double> >& q, int n_cells, double gamma) {

    // Calculate pressure from a given Q vector

    // Initialize pressure vector
    vector<double> p;
    p.resize(n_cells);

    // Calculate values
    for (int i=0; i<q[0].size(); i++) {
        p[i] = (gamma - 1)*(q.at(2).at(i) - q.at(1).at(i)*q.at(1).at(i)/(2*q.at(0).at(i)));
    }

    return p;

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
        this->iterate(this->q, this->time, this->gamma, this->dt, this->s, this->n_cells);
        this->write();
    }

}

void Flow::iterate(vector<vector<double> >& q, double& time, double& gamma, double& dt, double& s, int& n_cells) {

    // Copy q from current timestep into q_old array
    vector<vector<double> > q_old(q);

    // Find old f array
    vector<vector<double> > f_old = this->calculate_f_vector(q, n_cells, gamma);

    // Advance every grid point forward one timestep
    for (int i=1; i<(n_cells-1); i++) {
        q.at(0).at(i) = q_old.at(0).at(i) - (dt/(2*s))*(f_old.at(0).at(i+1) - f_old.at(0).at(i-1));
        q.at(1).at(i) = q_old.at(1).at(i) - (dt/(2*s))*(f_old.at(1).at(i+1) - f_old.at(1).at(i-1));
        q.at(2).at(i) = q_old.at(2).at(i) - (dt/(2*s))*(f_old.at(2).at(i+1) - f_old.at(2).at(i-1));
    }

    // Find pressure distribution
    this->p = this->calculate_pressure(q, n_cells, gamma);

    time = time + dt;

}

void Flow::write() {

    // Write data to solution.dat

    ofstream solution_file;
    solution_file.open("solution.dat", std::ios_base::app);
    solution_file << this->time << endl;
    for (int i = 0; i<this->n_cells; i++) {
        solution_file << this->grid[i] << " "
            << this->q.at(0).at(i) << " " << this->q.at(1).at(i) << " " << this->q.at(2).at(i) << " "
            << this->f.at(0).at(i) << " " << this->f.at(1).at(i) << " " << this->f.at(2).at(i) << " "
            << this->p.at(i) << endl;
    }
    solution_file << endl;

}

void Flow::output() {

    // Output results to stdout
    // TODO: Formatted output for a Python plotting script or Tecplot
    cout << "Final t = " << this->time << " s." << endl;
    cout << "Solution:" << endl;
    for (int i = 0; i < this->q[0].size(); i++) {
        cout << this->q.at(0).at(i) << ", " << this->q.at(1).at(i) << ", " << this->q.at(2).at(i) << "    ";
    }
    cout << endl;

}
