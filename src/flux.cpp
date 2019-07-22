#include <cmath>
#include <vector>
#include <flux.h>
using std::sqrt;
using std::vector;


// Calculate fluxes using Steger and Warming's method. Most of the math/theory
// behind this code is given in "Riemann Solvers and Numerical Methods for Fluid
// Dynamics" by E. F. Toro. The result is a vector with two elements, the first
// containing a vector of F_i-1/2, the second containing a vector of F_i+1/2.
vector<double> Flux::steger_warming(vector<double> q_left,
    vector<double> q_right, double gamma) {

    // Useful references and variables
    double &q1_l = q_left[0];
    double &q2_l = q_left[1];
    double &q3_l = q_left[2];
    double &q1_r = q_right[0];
    double &q2_r = q_right[1];
    double &q3_r = q_right[2];

    // Calculate speed of fluid, u, from conserved variables
    double u_l = q2_l/q1_l;
    double u_r = q2_r/q1_r;

    // Calculate speed of sound, a, from conserved variables
    // a = sqrt( gamma*R*T )
    // a = sqrt( gamma(gamma-1)(q3/q1 - (1/2)q2^2/q1) )
    double a_l = sqrt( gamma*(gamma-1)*(q3_l/q1_l - (1/2)*(q2_l*q2_l)/q1_l) );
    double a_r = sqrt( gamma*(gamma-1)*(q3_r/q1_r - (1/2)*(q2_r*q2_r)/q1_r) );

    // Calculate specific enthalpy, h, using a and conserved variables
    // h = (1/2)u^2 + (a^2)/(gamma-1)
    double h_l = (1/2)*(u_l*u_l) + (a_l*a_l)/(gamma-1);
    double h_r = (1/2)*(u_r*u_r) + (a_r*a_r)/(gamma-1);

    // Get eigenvalues
    vector<double> eigvalues_l;
    vector<double> eigvalues_r;
    eigvalues_l[0] = u_l - a_l;
    eigvalues_l[1] = u_l;
    eigvalues_l[2] = u_l + a_l;
    eigvalues_r[0] = u_r - a_r;
    eigvalues_r[1] = u_r;
    eigvalues_r[2] = u_r + a_r;

    // Split eigenvalues. Keep positive eigenvalues from the left, and negative
    // eigenvalues from the right.
    vector<double> eig_split_l;
    vector<double> eig_split_r;
    eig_split_l.reserve(3);
    eig_split_r.reserve(3);
    for (int i = 0; i < eigvalues_l.size(); i++) {
        if (eigvalues_l[i] > 0) {
            eig_split_l[i] = eigvalues_l[i];
        } else {
            eig_split_l[i] = 0;
        }
    }
    for (int i = 0; i < eigvalues_r.size(); i++) {
        if (eigvalues_r[i] < 0) {
            eig_split_r[i] = eigvalues_r[i];
        } else {
            eig_split_r[i] = 0;
        }
    }

    // Initialize and calculate fluxes
    vector<double> flux_l;
    vector<double> flux_r;
    flux_l.reserve(3);
    flux_r.reserve(3);
    double coefficient_l = q1_l/(2*gamma);
    double coefficient_r = q1_r/(2*gamma);
    // Calculate left fluxes
    flux_l[0] = coefficient_l * ( eig_split_l[0] +
        2*(gamma-1)*eig_split_l[1] + eig_split_l[2] );
    flux_l[1] = coefficient_l * ( (u_l-a_l)*eig_split_l[0] +
        2*(gamma-1)*u_l*eig_split_l[1] + (u_l+a_l)*eig_split_l[2] );
    flux_l[2] = coefficient_l * ( (h_l - u_l*a_l)*eig_split_l[0] +
        (gamma-1)*(u_l*u_l)*eig_split_l[1] + (h_l + u_l*a_l)*eig_split_l[2] );
    // Calculate right fluxes
    flux_r[0] = coefficient_r * ( eig_split_r[0] +
        2*(gamma-1)*eig_split_r[1] + eig_split_r[2] );
    flux_r[1] = coefficient_r * ( (u_r-a_r)*eig_split_r[0] +
        2*(gamma-1)*u_r*eig_split_r[1] + (u_r+a_r)*eig_split_r[2] );
    flux_r[2] = coefficient_r * ( (h_r - u_r*a_r)*eig_split_r[0] +
        (gamma-1)*(u_r*u_r)*eig_split_r[1] + (h_r + u_r*a_r)*eig_split_r[2] );

    // Add left and right fluxes to get resultant flux through face
    vector<double> total_flux;
    total_flux.reserve(3);
    for (int i = 0; i < total_flux.size(); i++) {
        total_flux[i] = flux_l[i] + flux_r[i];
    }

    // Return result
    return total_flux;

}