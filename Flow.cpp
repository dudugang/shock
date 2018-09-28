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
    this->n_cells = 50;
    this->dt = .001;
    this->time = 0;
    this->n_iter = 10;
    this->length = 10;
    // Left boundary
    this->boundary_conditions[0,0] = 1.225;    // rho
    this->boundary_conditions[1,0] = 10;       // u
    this->boundary_conditions[2,0] = 215400;   // e
    // Right boundary
    this->boundary_conditions[0,1] = 1.225;    // rho
    this->boundary_conditions[1,1] = 10;       // u
    this->boundary_conditions[2,1] = 215400;   // e

}


void Flow::initialize() {

    // Allocate grid and q arrays
    this->grid.resize(this->n_points);
    this->q[0].resize(this->n_points);
    this->q[1].resize(this->n_points);
    this->q[2].resize(this->n_points);

    // Generate grid
    for (int i=0; i<this->n_cells; i++) {
        this->grid[i] = (this->length*i)/(n_cells-2) - (this->length/2)/(n_cells-2);
    }
    // Calculate spacing
    this->s = this->length/(this->n_cells-2);

    // Generate initial conditions
    // TODO: Add user-input initial conditions to replace hard-coded values
    for (int i=1; i<(this->n_cells/3); i++) {
        this->q[0,i] = 1.225;
        this->q[1,i] = 10;
        this->q[2,i] = 215400;
    }
    for (int i=this->n_cells/3; i<(2*this->n_cells/3); i++) {
        this->q[0,i] = 2;
        this->q[1,i] = 40;
        this->q[2,i] = 300000;
    }
    for (int i=2*this->n_cells/3; i<this->n_cells; i++) {
        this->q[0,i] = 1.225;
        this->q[1,i] = 10;
        this->q[2,i] = 215400;
    }

    // Add boundary conditions
    this->q[0,0] = this->boundary_conditions[0,0]
    this->q[1,0] = this->boundary_conditions[1,0]
    this->q[2,0] = this->boundary_conditions[2,0]
    this->q[0,this->n_cells-1] = this->boundary_conditions[0,1]
    this->q[1,this->n_cells-1] = this->boundary_conditions[1,1]
    this->q[2,this->n_cells-1] = this->boundary_conditions[2,1]

}

void Flow::solve() {

    // Wipe old solution file
    remove("solution.dat");

    // Iterate until n_iter is reached
    for (int i=1; i<=this->n_iter; i++) {
        cout << "----- Iteration " << i << ", t = " << this->time << " s. " << endl;
        if (this->implicit == true) {
            this->iterateImplicit(this->u, this->time, this->alpha, this->dt, this->s);
            this->write();
        } else {
            this->iterateExplicit(this->u, this->time, this->alpha, this->dt, this->s);
        }
    }

}

void Flow::iterateExplicit(vector<double>& u, double& time, double alpha, double dt, double s) {

    // Copy u from current timestep into u_old array
    vector<double> u_old = u;

    // Advance every grid point forward one timestep
    for (int x=1; x<(u.size()-1); x++) {
        u[x] =
               u_old[x]
            + (alpha*dt)/(s*s)
            * (u_old[x+1] - 2.0*u_old[x] + u_old[x-1]);
    }

    time = time + dt;

}

void Flow::iterateImplicit(vector<double>& u, double& time, double alpha, double dt, double s) {

    // Create A matrix of linear system
    vector<vector<double>> A_matrix(u.size(), vector<double>(u.size(),0.0));
    A_matrix[0][0] = 1;
    A_matrix[u.size()-1][u.size()-1] = 1;
    for (int i=1; i<(u.size()-1); i++) {
        A_matrix[i][i] = 1 - (2*dt)/s;
        A_matrix[i][i+1] = dt/s;
        A_matrix[i][i-1] = dt/s;
    }

    // Create b matrix of linear system
    vector<double> b_matrix(u.size());
    for (int i=0; i<b_matrix.size(); i++) {
        b_matrix[i] = u[i];
    }

    // Solve linear system
    vector<double> solution = this->solveLinearSystem(A_matrix, b_matrix);

    // Update u vector
    for (int i=0; i<(u.size()); i++) {
        u[i] = solution[i];
    }

    time = time + dt;

}

vector<double> Flow::solveLinearSystem(vector<vector<double>>& A_matrix, vector<double>& b_matrix) {

    // Solve linear system using the Jacobi method
    // TODO: Create separate class for solving linear systems

    // Create remainder and diagonal inverse matrices
    vector<vector<double>> R_matrix = A_matrix;
    vector<vector<double>> D_inverse_matrix(A_matrix.size(), vector<double>(A_matrix.size()));
    for (int i=0; i<A_matrix.size(); i++) {
        R_matrix[i][i] = 0;
        D_inverse_matrix[i][i] = 1/A_matrix[i][i];
    }

    // Iterate over system
    // TODO: Check for error to know when to stop iterating
    vector<double> x(A_matrix.size(),0.0);
    for (int i=0; i<100; i++) {
        vector<double> R_times_x = this->multiplyMatrices(R_matrix, x);
        vector<double> b_minus_Rx = this->subtractMatrices(b_matrix, R_times_x);
        x = this->multiplyMatrices(D_inverse_matrix, b_minus_Rx);
    }

    return x;

}

vector<double> Flow::multiplyMatrices(vector<vector<double>>& left, vector<double>& right) {

    // NOTE: This function can only multiply a 2D matrix with a 1D matrix

    vector<double> result(left.size());

    for (int i=0; i<left.size(); i++) {
        for (int j=0; j<right.size(); j++) {
            result[i] = result[i] + left[i][j]*right[j];
        }
    }

    return result;

}

vector<double> Flow::subtractMatrices(vector<double>& left, vector<double>& right) {

    // NOTE: This function can only subtract two 1D matrices

    vector<double> result(left.size());

    for (int i=0; i<left.size(); i++) {
        result[i] = left[i] - right[i];
    }

    return result;

}

void Flow::write() {

    // Write data to solution.dat
    ofstream solution_file;
    solution_file.open("solution.dat", std::ios_base::app);
    solution_file << this->time << endl;
    for (int i = 0; i<this->u.size(); i++) {
        solution_file << i*this->s << "    " << this->u[i] << endl;
    }
    solution_file << endl;

}

void Flow::output() {

    // Output results to stdout
    // TODO: Formatted output for a Python plotting script or Tecplot
    cout << "Final t = " << this->time << " s." << endl;
    cout << "Solution:" << endl;
    for (int i = 0; i < this->u.size(); i++) {
        cout << this->u[i] << " ";
    }
    cout << endl;

}
