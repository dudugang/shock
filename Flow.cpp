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
    this->n_cells = 100;
    this->cfl = .1;
    this->time = 0;
    this->n_iter = 2000;
    this->length = 1;
    // Left boundary
    this->boundary_conditions.at(0).at(0) = 1;    // rho
    this->boundary_conditions.at(1).at(0) = 1;    // rho u
    this->boundary_conditions.at(2).at(0) = 120;  // rho e
    // Right boundary
    this->boundary_conditions.at(0).at(1) = 1;    // rho
    this->boundary_conditions.at(1).at(1) = 1;    // rho u
    this->boundary_conditions.at(2).at(1) = 100;  // rho e
    this->gamma = 1.4;
}


void Flow::initialize() {
    // Allocate grid, q, q_vertex, f, and p vectors
    this->grid.resize(this->n_cells);
    this->q.resize(3);
    this->q[0].resize(this->n_cells);
    this->q[1].resize(this->n_cells);
    this->q[2].resize(this->n_cells);
    this->q_vertex.resize(3);
    this->q_vertex[0].resize(this->n_cells-1);
    this->q_vertex[1].resize(this->n_cells-1);
    this->q_vertex[2].resize(this->n_cells-1);
    this->f.resize(3);
    this->f[0].resize(this->n_cells);
    this->f[1].resize(this->n_cells);
    this->f[2].resize(this->n_cells);
    this->p.resize(this->n_cells);

    // Generate grid
    for (int i=0; i<this->n_cells; i++) {
        this->grid.at(i) = (this->length*i)/(n_cells-2) - (this->length/2)/(n_cells-2);
    }

    // Calculate spacing
    this->s = this->length/(this->n_cells-2);

    // Generate initial conditions
    // TODO: Add user-input initial conditions to replace hard-coded values

    // Constant
    for (int i = 1; i<(this->n_cells)-1; i++) {
        this->q.at(0).at(i) = 1;
        this->q.at(1).at(i) = 1;
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
    this->q.at(0).at(0) = this->boundary_conditions.at(0).at(0);
    this->q.at(1).at(0) = this->boundary_conditions.at(1).at(0);
    this->q.at(2).at(0) = this->boundary_conditions.at(2).at(0);
    this->q.at(0).at(this->n_cells-1) = this->boundary_conditions.at(0).at(1);
    this->q.at(1).at(this->n_cells-1) = this->boundary_conditions.at(1).at(1);
    this->q.at(2).at(this->n_cells-1) = this->boundary_conditions.at(2).at(1);

    // Calculate cell wall Q's (average neighbor cells)
    for (int i=0; i<this->q_vertex[0].size(); i++) {
        q_vertex.at(0).at(i) = (q.at(0).at(i) + q.at(0).at(i+1))/2;
        q_vertex.at(1).at(i) = (q.at(1).at(i) + q.at(1).at(i+1))/2;
        q_vertex.at(2).at(i) = (q.at(2).at(i) + q.at(2).at(i+1))/2;
    }

    // Calculate timestep
    this->calculate_dt(this->q, this->gamma, this->cfl, this->s);
}

vector<vector<double> > Flow::calculate_f_vector(vector<vector<double> >& q, int n_cells, double gamma) {
    // Calculate F vector from a given Q vector

    // Initialize F vector
    vector<vector<double> > f;
    f.resize(3);
    f[0].resize(q[0].size());
    f[1].resize(q[0].size());
    f[2].resize(q[0].size());

    // Calculate values
    for (int i=0; i<q[0].size(); i++) {
        f.at(0).at(i) = q.at(1).at(i);
        f.at(1).at(i) = q.at(1).at(i)*q.at(1).at(i)/q.at(0).at(i) + (gamma - 1)*(q.at(2).at(i) - q.at(1).at(i)*q.at(1).at(i)/(2*q.at(0).at(i)));
        f.at(2).at(i) = q.at(1).at(i)*q.at(2).at(i)/q.at(0).at(i) + (gamma - 1)*(q.at(1).at(i)/q.at(0).at(i))*(q.at(2).at(i) - q.at(1).at(i)*q.at(1).at(i)/(2*q.at(0).at(i)));
    }

    this->f = f;
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
        this->iterate(this->q, this->q_vertex, this->time, this->gamma, this->dt, this->s, this->n_cells);
        this->write();
    }
}


// // // // //
// This function runs one iteration of the flow solver. This involves
// calculating wall fluxes by solving the Riemann problem, calculating the
// timestep given the CFL, and advancing the solution forward to the next
// timestep.
void Flow::iterate(
        vector<vector<double> >& q, vector<vector<double> >& q_vector,
        double& time, double& gamma, double& cfl, double& s, int& n_cells) {

    // Copy q from current timestep into q_old array
    vector<vector<double> > q_old(q);

    // Write this function
    // f_old = this->calculate_riemann_fluxes

    // Find old f array
    vector<vector<double> > f_old = this->calculate_f_vector(q_vertex, n_cells, gamma);

    // Calculate timestep
    dt = this->calculate_dt(this->q, this->gamma, this->cfl, this->s);

    cout << "F1: " << f_old.at(0).at(1) << endl;
    // Advance every grid point forward one timestep
    for (int i=1; i<(n_cells-1); i++) {
        q.at(0).at(i) = q_old.at(0).at(i) - (dt/(s))*(f_old.at(0).at(i) - f_old.at(0).at(i-1));
        q.at(1).at(i) = q_old.at(1).at(i) - (dt/(s))*(f_old.at(1).at(i) - f_old.at(1).at(i-1));
        q.at(2).at(i) = q_old.at(2).at(i) - (dt/(s))*(f_old.at(2).at(i) - f_old.at(2).at(i-1));
    }

    // Update vertex grid
    // TODO: Make this (as well as initialization) into a function
    for (int i=0; i<this->q_vertex[0].size(); i++) {
        q_vertex.at(0).at(i) = (q.at(0).at(i) + q.at(0).at(i+1))/2;
        q_vertex.at(1).at(i) = (q.at(1).at(i) + q.at(1).at(i+1))/2;
        q_vertex.at(2).at(i) = (q.at(2).at(i) + q.at(2).at(i+1))/2;
    }

    // Find pressure distribution
    this->p = this->calculate_pressure(q, gamma);

    time = time + dt;
}

// // // // //
// This function calculates p_star, which is the pressure on either side of the
// contact surface when solving the general Riemann problem. The equations for
// this are given in Knight, and the solution is done using Newton's Method.
double Flow::calculate_pstar(
        double p1, double p4, double u1, double u4, double a1, double a4,
        double gamma) {

    // TODO: Make a smarter initial guess (Knight gives a good one)
    double pstar = p1;
    // TODO: Have iteration stop by an error tolerance as opposed to a
    // hardcoded number of iterations
    int newton_iters = 20;

    double function;
    double function_prime;
    for (int i = 0; i < newton_iters; i++) {
        function = a1*f_pstar(pstar, p1, gamma) + a4*f_pstar(pstar, p4, gamma) - u1 + u4;
        function_prime = a1*f_prime_pstar(pstar, p1, gamma) + a4*f_prime_pstar(pstar, p4, gamma);
        pstar = pstar - function/function_prime;
    }

    return pstar;
}

// // // // //
// This function calculates the function f(p_star, p) defined in Knight, page 33,
// eq. 2.99. This is used in the calculation of p_star when solving the general
// Reimann problem at every cell interface.
double Flow::f_pstar(
        double ps, double p, double gamma) {

    if (ps >= p) {
        return ( (ps/p - 1)/(gamma*pow((gamma - 1)/(2*gamma) + (ps*(gamma + 1))/(2*gamma*p), 1/2)) );
    } else {
        return ( (2*(pow((ps/p), (gamma - 1)/(2*gamma)) - 1))/(gamma - 1) );
    }
}

// // // // //
// This function calculates the derivative of the function f(p_star, p) defined
// in Knight, page 33, eq. 2.99. This is used in the calculation of p_star for
// iterations of Newton's method.
double Flow::f_prime_pstar(double ps, double p, double gamma) {

    if (ps >= p) {
        return ( (pow(2, 1/2)*(ps - p + 3*gamma*p + gamma*ps))/(2*pow(gamma, 2)*pow(p, 2)*pow((- p + ps + gamma*p + gamma*ps)/(gamma*p), 3/2)) );
    } else {
        return 1/(gamma*p*pow(ps/p, (gamma + 1)/(2*gamma)));
    }
}

void Flow::write() {

    // Write data to solution.dat

    ofstream solution_file;
    solution_file.open("solution.dat", std::ios_base::app);
    solution_file << this->time << endl;
    for (int i = 0; i<this->n_cells; i++) {
        solution_file << this->grid[i] << " "
            << this->q.at(0).at(i) << " " << this->q.at(1).at(i) << " " << this->q.at(2).at(i) << " "
            << this->p.at(i) << endl;
            //<< this->f.at(0).at(i) << " " << this->f.at(1).at(i) << " " << this->f.at(2).at(i) << " "
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
