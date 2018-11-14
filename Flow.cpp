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
    this->n_cells = 500;
    this->cfl = .2;
    this->time = 0;
    this->n_iter = 1000;
    this->length = 1;
    // Left boundary
    this->boundary_conditions.at(0).at(0) = 1;    // rho
    this->boundary_conditions.at(1).at(0) = 1;    // rho u
    this->boundary_conditions.at(2).at(0) = 100;  // rho e
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
    this->q_vertex[0].resize(this->n_cells+1);
    this->q_vertex[1].resize(this->n_cells+1);
    this->q_vertex[2].resize(this->n_cells+1);
    this->f.resize(3);
    this->f[0].resize(this->n_cells+1);
    this->f[1].resize(this->n_cells+1);
    this->f[2].resize(this->n_cells+1);
    this->p.resize(this->n_cells);

    // Generate grid
    for (int i=0; i<this->n_cells; i++) {
        this->grid.at(i) = (this->length*i)/(n_cells-2) - (this->length/2)/(n_cells-2);
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
        this->q.at(1).at(i) = 1;
        this->q.at(2).at(i) = 150;
    }
    for (int i = this->n_cells/2; i < this->n_cells; i++) {
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

    // Update q_vertex by solving the Reimann problem at every cell interface
    for (int i = 1; i < n_cells; i++) {
        // This is eq. 2.4 from Knight
        double p1 = (gamma-1)*( q.at(2).at(i-1) - pow(q.at(1).at(i-1), 2)/(2*q.at(0).at(i-1)) );
        double p4 = (gamma-1)*( q.at(2).at(i) - pow(q.at(1).at(i), 2)/(2*q.at(0).at(i)) );
        // This is because q = (rho, rho*u, rho*e)
        double u1 = q.at(1).at(i-1) / q.at(0).at(i-1);
        double u4 = q.at(1).at(i) / q.at(0).at(i);
        // This comes from a = sqrt(gamma*p/rho)
        double a1 = pow(gamma*p1/q.at(0).at(i-1), .5);
        double a4 = pow(gamma*p4/q.at(0).at(i), .5);

        // Calculate pstar for the Riemann problem
        double pstar = this->calculate_pstar(
                p1, p4, u1, u4, a1, a4, gamma);

        // Use this ptar to solve the Riemann problem
        vector<double> rho_u_p = this->solve_riemann_problem(
                p1, p4, u1, u4, a1, a4, gamma, pstar);

        // Convert rho, u, p to Q vertex values
        q_vertex.at(0).at(i) = rho_u_p[0];
        q_vertex.at(1).at(i) = rho_u_p[0] * rho_u_p[1];
        q_vertex.at(2).at(i) = rho_u_p[2] / (gamma-1);
    }
    // TEMPORARY: Reflective boundary
    q_vertex.at(0).at(0) = q_vertex.at(0).at(1);
    q_vertex.at(1).at(0) = q_vertex.at(1).at(1);
    q_vertex.at(2).at(0) = q_vertex.at(2).at(1);
    q_vertex.at(0).at(n_cells) = q_vertex.at(0).at(n_cells-1);
    q_vertex.at(1).at(n_cells) = q_vertex.at(1).at(n_cells-1);
    q_vertex.at(2).at(n_cells) = q_vertex.at(2).at(n_cells-1);

    // Find current f array
    for (int i = 0; i < (n_cells + 1); i++) {
        vector<double> q_interface;
        q_interface.resize(3);
        vector<double> f_interface;
        f_interface.resize(3);
        q_interface[0] = q_vertex.at(0).at(i);
        q_interface[1] = q_vertex.at(1).at(i);
        q_interface[2] = q_vertex.at(2).at(i);
        f_interface = this->calculate_f_vector(q_interface, gamma);
        f.at(0).at(i) = f_interface[0];
        f.at(1).at(i) = f_interface[1];
        f.at(2).at(i) = f_interface[2];
    }

    // Calculate timestep
    dt = this->calculate_dt(this->q, this->gamma, this->cfl, this->s);

    // Advance every grid point forward one timestep
    for (int i = 0; i < n_cells; i++) {
        q.at(0).at(i) = q.at(0).at(i) - (dt/(s))*(f.at(0).at(i+1) - f.at(0).at(i));
        q.at(1).at(i) = q.at(1).at(i) - (dt/(s))*(f.at(1).at(i+1) - f.at(1).at(i));
        q.at(2).at(i) = q.at(2).at(i) - (dt/(s))*(f.at(2).at(i+1) - f.at(2).at(i));
    }

    // Find pressure distribution
    this->p = this->calculate_pressure(q, gamma);

    time = time + dt;
}

// // // // //
// This function solves the General Riemann Problem, as defined by Knight in
// chapter 2.9. The individual cases are written in seperate functions, so this
// function just picks which case to run given the inputs.
vector<double> Flow::solve_riemann_problem(
        double p1, double p4, double u1, double u4, double a1, double a4,
        double gamma, double pstar) {

    // Initialize vector of results
    vector<double> results;

    // The following conditionals are defined by Knight in chapter 2.9
    if ((p1 <= pstar) && (p4 <= pstar)) {
        results = this->case_1_riemann(p1, p4, u1, u4, a1, a4, gamma, pstar);
    } else if ((p1 <= pstar) && (p4 >= pstar)) {
        results = this->case_2_riemann(p1, p4, u1, u4, a1, a4, gamma, pstar);
    } else if ((p1 >= pstar) && (p4 <= pstar)) {
        results = this->case_3_riemann(p1, p4, u1, u4, a1, a4, gamma, pstar);
    } else if ((p1 >= pstar) && (p4 >= pstar)) {
        results = this->case_4_riemann(p1, p4, u1, u4, a1, a4, gamma, pstar);
    } else {
        cout << "Riemann problem failed to solve." << endl;
        cout << "pstar: " << pstar << endl;
        cout << "p1, u1, a1: " << p1 << ", " << u1 << ", " << a1 << endl;
        cout << "p4, u4, a4: " << p4 << ", " << u4 << ", " << a4 << endl;
    }

    return results;
}

// // // // //
// This function calculates the flow variables at the cell interface given by
// Case 1 of the general Riemann problem, as defined by Knight. First, the
// velocities of each characteristic are calculated. This information is then
// used to pick the correct region of the x-t diagram from which to extract flow
// information.
vector<double> Flow::case_1_riemann(
        double p1, double p4, double u1, double u4, double a1, double a4,
        double gamma, double pstar) {

    // Initialize vector of results, which will contain the desired (rho, u, p)
    vector<double> results;
    results.resize(3);

    // Find rho1 and rho4 using a = sqrt(gamma*p/rho)
    double rho1 = gamma*p1/(pow(a1,2));
    double rho4 = gamma*p4/(pow(a4,2));

    // Calculate left shock velocity
    double c_shock_left = u1 - a1*pow( ((gamma+1)/(2*gamma)) * (pstar/p1 - 1) + 1, .5);

    // Calculate right shock velocity
    double c_shock_right = u4 + a4*pow( ((gamma+1)/(2*gamma)) * (pstar/p4 - 1) + 1, .5);

    // Calculate contact surface velocity
    double c_contact_surface = u1 - (a1/gamma)*(pstar/p1 - 1)*pow( ((gamma+1)/(2*gamma))*(pstar/p1) + (gamma-1)/(2*gamma), -.5);

    // If the leftmost characteristic is going to the right, then the flow
    // properties at the cell wall are equal to state 1. The opposite is also
    // true; if the rightmost characteristic is going to the left, then the flow
    // properties at the cell wall are equal to state 4.
    if (c_shock_left >= 0) {
        results[0] = rho1;
        results[1] = u1;
        results[2] = p1;
    } else if (c_shock_right <= 0) {
        results[0] = rho4;
        results[1] = u4;
        results[2] = p4;
    }
    // If we've reached this point, then we know that the leftmost
    // characteristic goes to the left, and the rightmost characteristic goes to
    // the right. Thus, the desired state lies between. If the contact surface
    // goes to the right, then the wall flow properties are state 2. If the
    // contact surface goes to the left, then the wall flow properties are state
    // 3. The following equations are defined in Knight.
      else if (c_contact_surface >= 0) {
        results[0] = rho1*( (gamma-1) + (gamma+1)*(pstar/p1) )/( (gamma+1) + (gamma-1)*(pstar/p1) );
        results[1] = c_contact_surface;
        results[2] = pstar;
    } else if (c_contact_surface <= 0) {
        results[0] = rho4*( (gamma-1) + (gamma+1)*(pstar/p4) )/( (gamma+1) + (gamma-1)*(pstar/p4) );
        results[1] = c_contact_surface;
        results[2] = pstar;
    } else {
        cout << "Case 1 of Reimann problem failed to solve." << endl;
    }

    return results;
}

// // // // //
// This function calculates the flow variables at the cell interface given by
// Case 2 of the general Riemann problem, as defined by Knight. First, the
// velocities of each characteristic are calculated. This information is then
// used to pick the correct region of the x-t diagram from which to extract flow
// information.
vector<double> Flow::case_2_riemann(
        double p1, double p4, double u1, double u4, double a1, double a4,
        double gamma, double pstar) {

    // Initialize vector of results, which will contain the desired (rho, u, p)
    vector<double> results;
    results.resize(3);

    // Find rho1 and rho4 using a = sqrt(gamma*p/rho)
    double rho1 = gamma*p1/(pow(a1,2));
    double rho4 = gamma*p4/(pow(a4,2));

    // Calculate left shock velocity
    double c_shock_left = u1 - a1*pow( ((gamma+1)/(2*gamma)) * (pstar/p1 - 1) + 1, .5);

    // Calculate velocity of left and right side of right expansion
    double c_expansion_right_left = u4 + a4*( ((gamma+1)/(gamma-1))*pow((pstar/p4), (gamma-1)/(2*gamma)) - 2/(gamma-1) );
    double c_expansion_right_right = u4 + a4;

    // Calculate contact surface velocity
    double c_contact_surface = u1 - (a1/gamma)*(pstar/p1 - 1)*pow( ((gamma+1)/(2*gamma))*(pstar/p1) + (gamma-1)/(2*gamma), -.5);

    // If the leftmost characteristic is going to the right, then the flow
    // properties at the cell wall are equal to state 1. The opposite is also
    // true; if the rightmost characteristic is going to the left, then the flow
    // properties at the cell wall are equal to state 4.
    if (c_shock_left >= 0) {
        results[0] = rho1;
        results[1] = u1;
        results[2] = p1;
    } else if (c_expansion_right_right <= 0) {
        results[0] = rho4;
        results[1] = u4;
        results[2] = p4;
    }
    // If we've reached this point, then we know that the leftmost
    // characteristic goes to the left, and the rightmost characteristic goes to
    // the right. Thus, the desired state lies between. If the contact surface
    // goes to the right, then the wall flow properties are state 2. If the
    // contact surface goes to the left and the left side of the expansion goes
    // to the right, then the wall flow properties are state 3. If the left side
    // of the expansion goes to the left and the right side of the expansion
    // goes to the right, then the flow properties are those within the fan,
    // with x/t = 0. The following equations are defined in Knight.
      else if (c_contact_surface >= 0) {
        results[0] = rho1*( (gamma-1) + (gamma+1)*(pstar/p1) )/( (gamma+1) + (gamma-1)*(pstar/p1) );
        results[1] = c_contact_surface;
        results[2] = pstar;
    } else if ((c_contact_surface <= 0) && (c_expansion_right_left >= 0)) {
        results[0] = rho4*pow(pstar/p4, 1/gamma);
        results[1] = c_contact_surface;
        results[2] = pstar;
    } else if ((c_expansion_right_left <= 0) && (c_expansion_right_right >= 0)) {
        results[1] = (2/(gamma+1))*( ((gamma-1)/2)*u4 - a4 );
        results[2] = p4*pow( ((gamma-1)/(gamma+1))*(1/a4)*(-u4) + 2/(gamma+1), (2*gamma)/(gamma-1) );
        results[0] = rho4*pow(results[2]/p4, 1/gamma);
    } else {
        cout << "Case 2 of Reimann problem failed to solve." << endl;
    }

    return results;
}

// // // // //
// This function calculates the flow variables at the cell interface given by
// Case 3 of the general Riemann problem, as defined by Knight. First, the
// velocities of each characteristic are calculated. This information is then
// used to pick the correct region of the x-t diagram from which to extract flow
// information.
vector<double> Flow::case_3_riemann(
        double p1, double p4, double u1, double u4, double a1, double a4,
        double gamma, double pstar) {

    // Initialize vector of results, which will contain the desired (rho, u, p)
    vector<double> results;
    results.resize(3);

    // Find rho1 and rho4 using a = sqrt(gamma*p/rho)
    double rho1 = gamma*p1/(pow(a1,2));
    double rho4 = gamma*p4/(pow(a4,2));

    // Calculate velocity of left and right side of left expansion
    double c_expansion_left_left = u1 - a1;
    double c_expansion_left_right = u1 + a1*( -((gamma+1)/(gamma-1))*pow((pstar/p1), (gamma-1)/(2*gamma)) + 2/(gamma-1) );

    // Calculate right shock velocity
    double c_shock_right = u4 + a4*pow( ((gamma+1)/(2*gamma)) * (pstar/p4 - 1) + 1, .5);

    // Calculate contact surface velocity
    double c_contact_surface = u4 + (a4/gamma)*(pstar/p4 - 1)*pow( ((gamma+1)/(2*gamma))*(pstar/p4) + (gamma-1)/(2*gamma), -.5);

    // If the leftmost characteristic is going to the right, then the flow
    // properties at the cell wall are equal to state 1. The opposite is also
    // true; if the rightmost characteristic is going to the left, then the flow
    // properties at the cell wall are equal to state 4.
    if (c_expansion_left_left >= 0) {
        results[0] = rho1;
        results[1] = u1;
        results[2] = p1;
    } else if (c_shock_right <= 0) {
        results[0] = rho4;
        results[1] = u4;
        results[2] = p4;
    }
    // If we've reached this point, then we know that the leftmost
    // characteristic goes to the left, and the rightmost characteristic goes to
    // the right. Thus, the desired state lies between. If the contact surface
    // goes to the left, then the wall flow properties are state 3. If the
    // contact surface goes to the right and the right side of the expansion goes
    // to the left, then the wall flow properties are state 2. If the left side
    // of the expansion goes to the left and the right side of the expansion
    // goes to the right, then the flow properties are those within the fan,
    // with x/t = 0. The following equations are defined in Knight.
      else if (c_contact_surface <= 0) {
        results[0] = rho4*( (gamma-1) + (gamma+1)*(pstar/p4) )/( (gamma+1) + (gamma-1)*(pstar/p4) );
        results[1] = c_contact_surface;
        results[2] = pstar;
    } else if ((c_contact_surface >= 0) && (c_expansion_left_right <= 0)) {
        results[0] = rho1*pow(pstar/p1, 1/gamma);
        results[1] = c_contact_surface;
        results[2] = pstar;
    } else if ((c_expansion_left_left <= 0) && (c_expansion_left_right >= 0)) {
        results[1] = (2/(gamma+1))*( ((gamma-1)/2)*u1 + a1 );
        results[2] = p1*pow( ((gamma-1)/(gamma+1))*(1/a1)*(u1) + 2/(gamma+1), (2*gamma)/(gamma-1) );
        results[0] = rho1*pow(results[2]/p1, 1/gamma);
    } else {
        cout << "Case 3 of Reimann problem failed to solve." << endl;
    }

    return results;
}

// // // // //
// This function calculates the flow variables at the cell interface given by
// Case 4 of the general Riemann problem, as defined by Knight. First, the
// velocities of each characteristic are calculated. This information is then
// used to pick the correct region of the x-t diagram from which to extract flow
// information.
vector<double> Flow::case_4_riemann(
        double p1, double p4, double u1, double u4, double a1, double a4,
        double gamma, double pstar) {

    // Initialize vector of results, which will contain the desired (rho, u, p)
    vector<double> results;
    results.resize(3);

    // Find rho1 and rho4 using a = sqrt(gamma*p/rho)
    double rho1 = gamma*p1/(pow(a1,2));
    double rho4 = gamma*p4/(pow(a4,2));

    // Calculate velocity of left and right side of left expansion
    double c_expansion_left_left = u1 - a1;
    double c_expansion_left_right = u1 + a1*( -((gamma+1)/(gamma-1))*pow((pstar/p1), (gamma-1)/(2*gamma)) + 2/(gamma-1) );

    // Calculate velocity of left and right side of right expansion
    double c_expansion_right_left = u4 + a4*( ((gamma+1)/(gamma-1))*pow((pstar/p4), (gamma-1)/(2*gamma)) - 2/(gamma-1) );
    double c_expansion_right_right = u4 + a4;

    // Calculate contact surface velocity
    double c_contact_surface = u1 + (2*a1/(gamma-1))*( 1 - pow(pstar/p1, (gamma-1)/(2*gamma)) );

    // If the leftmost characteristic is going to the right, then the flow
    // properties at the cell wall are equal to state 1. The opposite is also
    // true; if the rightmost characteristic is going to the left, then the flow
    // properties at the cell wall are equal to state 4.
    if (c_expansion_left_left >= 0) {
        results[0] = rho1;
        results[1] = u1;
        results[2] = p1;
    } else if (c_expansion_right_right <= 0) {
        results[0] = rho4;
        results[1] = u4;
        results[2] = p4;
    }
    // If we've reached this point, then we know that the leftmost
    // characteristic goes to the left, and the rightmost characteristic goes to
    // the right. Thus, the desired state lies between. If the contact surface
    // goes to the left and the left side of the right expansion goes to the
    // right, then the wall flow properties are state 3. If the contact surface
    // goes to the right and the right side of the left expansion goes to the
    // left, then the wall flow properties are state 2. If the left side of the
    // left expansion goes to the left and the right side of the left expansion
    // goes to the right, then the flow properties are those within the left
    // fan, with x/t = 0. if the left side of the right expansion goes to the
    // left and the right side of the right expansion goes to the right, then
    // the flow properties are those within the right fan, with x/t = 0. The
    // following equations are defined in Knight.
      else if ((c_contact_surface <= 0) && (c_expansion_right_left >= 0)) {
        results[0] = rho4*pow(pstar/p4, 1/gamma);
        results[1] = c_contact_surface;
        results[2] = pstar;
    } else if ((c_contact_surface >= 0) && (c_expansion_left_right <= 0)) {
        results[0] = rho1*pow(pstar/p1, 1/gamma);
        results[1] = c_contact_surface;
        results[2] = pstar;
    } else if ((c_expansion_left_left <= 0) && (c_expansion_left_right >= 0)) {
        results[1] = (2/(gamma+1))*( ((gamma-1)/2)*u1 + a1 );
        results[2] = p1*pow( ((gamma-1)/(gamma+1))*(1/a1)*(u1) + 2/(gamma+1), (2*gamma)/(gamma-1) );
        results[0] = rho1*pow(results[2]/p1, 1/gamma);
    } else if ((c_expansion_right_left <= 0) && (c_expansion_right_right >= 0)) {
        results[1] = (2/(gamma+1))*( ((gamma-1)/2)*u4 - a4 );
        results[2] = p4*pow( ((gamma-1)/(gamma+1))*(1/a4)*(-u4) + 2/(gamma+1), (2*gamma)/(gamma-1) );
        results[0] = rho4*pow(results[2]/p4, 1/gamma);
    } else {
        cout << "Case 4 of Reimann problem failed to solve." << endl;
    }

    return results;
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
