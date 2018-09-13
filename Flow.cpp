#include <array>
#include <iostream>
#include <vector>
#include "Flow.h"
using std::array;
using std::cout;
using std::endl;
using std::vector;

Flow::Flow() {

    // Generate parameters and boundary conditions
    // TODO: Add inputfile so that inputs are not hard-coded
    // TODO: Add adiabatic boundary condition
    this->n_points = 5;
    this->dt = 100;
    this->time = 0;
    this->n_iter = 1000;
    this->length = 10;
    this->alpha = 1;
    this->boundary_u[0] = 500;
    this->boundary_u[1] = 100;
    this->implicit = true;

}


void Flow::initialize() {

    // Allocate grid and u arrays
    this->grid.resize(this->n_points);
    this->u.resize(this->n_points);

    // Generate grid
    for (int i=0; i<this->n_points; i++) {
        this->grid[i] = ((i)/(this->n_points - 1)) * this->length;
    }
    // Calculate spacing
    this->s = this->length/(this->n_points-1);

    // Generate initial conditions
    // TODO: Add user-input initial conditions to replace hard-coded values
    for (int i=1; i<(this->n_points/2); i++) {
        this->u[i] = 500;
    }
    for (int i=this->n_points/2; i<this->u.size()-1; i++) {
        this->u[i] = 100;
    }

    // Add boundary conditions
    this->u[0] = this->boundary_u[0];
    this->u[this->u.size()-1] = this->boundary_u[1];

}

void Flow::solve() {

    // Iterate until n_iter is reached
    for (int i=1; i<=this->n_iter; i++) {
        cout << "----- Iteration " << i << ", t = " << this->time << " s. " << endl;
        if (this->implicit == true) {
            this->iterateImplicit(this->u, this->time, this->alpha, this->dt, this->s);
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
    vector<vector<double>> A_matrix(u.size()-2, vector<double>(u.size()-2,0.0));
    for (int i=0; i<(u.size()-2); i++) {
        A_matrix[i][i] = ((s*s)/(alpha*dt)) + 2;
    }
    for (int i=0; i<(u.size()-3); i++) {
        A_matrix[i][i+1] = -1;
        A_matrix[i+1][i] = -1;
    }

    // Create b matrix of linear system
    vector<double> b_matrix(u.size()-2);
    b_matrix[0] = ((s*s)/(alpha*dt))*u[1] + u[0];
    b_matrix[b_matrix.size()-1] = ((s*s)/(alpha*dt))*u[u.size()-2] + u[u.size()-1];
    for (int i=1; i<(b_matrix.size()-1); i++) {
        b_matrix[i] = (s*s)/(alpha*dt)*u[i+1];
    }

    // Solve linear system
    vector<double> solution = this->solveLinearSystem(A_matrix, b_matrix);

    // Update u vector
    for (int i=1; i<(u.size()-1); i++) {
        u[i] = solution[i-1];
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
