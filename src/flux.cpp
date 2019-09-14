#include <flux.h>
#include <state.h>


// Constructor
Flux::Flux() {

    // Resize vectors
    eigvalues_l.resize(4);
    eigvalues_r.resize(4);
    eig_plus.resize(4);
    eig_minus.resize(4);
    flux_plus.resize(4);
    flux_minus.resize(4);
    total_flux.resize(4);

}


// Calculate fluxes using Steger and Warming's method. Most of the math/theory
// behind this code is given in "Riemann Solvers and Numerical Methods for Fluid
// Dynamics" by E. F. Toro. The result is a vector with two elements, the first
// containing a vector of F_i-1/2, the second containing a vector of F_i+1/2.
vector<double> Flux::steger_warming(vector<double> q_left,
    vector<double> q_right, double gamma) {

    // Convert conserved variables to physical variables
    vector<double> physical_l = State::conserved_to_physical(q_left, gamma);
    vector<double> physical_r = State::conserved_to_physical(q_right, gamma);

    // Useful references
    double rho_l = physical_l[0];
    double u_l   = physical_l[1];
    double v_l   = physical_l[2];
    double p_l   = physical_l[3];
    double rho_r = physical_r[0];
    double u_r   = physical_r[1];
    double v_r   = physical_r[2];
    double p_r   = physical_r[3];

    // Calculate speed of sound, a, from conserved variables
    // a = sqrt( gamma*p/rho )
    double a_l = sqrt( gamma*p_l/rho_l );
    double a_r = sqrt( gamma*p_r/rho_r );

    // Calculate specific enthalpy, h, using a and conserved variables
    // h = (1/2)V^2 + (a^2)/(gamma-1)
    double h_l = (.5)*(u_l*u_l + v_l*v_l) + (a_l*a_l)/(gamma-1);
    double h_r = (.5)*(u_r*u_r + v_r*v_r) + (a_r*a_r)/(gamma-1);

    // Get eigenvalues
    eigvalues_l[0] = u_l - a_l;
    eigvalues_l[1] = u_l;
    eigvalues_l[2] = u_l;
    eigvalues_l[3] = u_l + a_l;
    eigvalues_r[0] = u_r - a_r;
    eigvalues_r[1] = u_r;
    eigvalues_r[2] = u_r;
    eigvalues_r[3] = u_r + a_r;

    // Split eigenvalues. Keep positive eigenvalues from the left, and negative
    // eigenvalues from the right.
    for (unsigned int i = 0; i < eigvalues_l.size(); i++) {
        if (eigvalues_l[i] > 0) {
            eig_plus[i] = eigvalues_l[i];
        } else {
            eig_plus[i] = 0;
        }
    }
    for (unsigned int i = 0; i < eigvalues_r.size(); i++) {
        if (eigvalues_r[i] < 0) {
            eig_minus[i] = eigvalues_r[i];
        } else {
            eig_minus[i] = 0;
        }
    }

    // Initialize and calculate fluxes
    double coefficient_l = rho_l/(2*gamma);
    double coefficient_r = rho_r/(2*gamma);
    // Calculate left fluxes
    flux_plus[0] = coefficient_l * ( eig_plus[0]
        + 2*(gamma-1)*eig_plus[1] + eig_plus[3] );
    flux_plus[1] = coefficient_l * ( (u_l-a_l)*eig_plus[0]
        + 2*(gamma-1)*u_l*eig_plus[1] + (u_l+a_l)*eig_plus[3] );
    flux_plus[2] = coefficient_l * ( v_l*eig_plus[0]
        + 2*(gamma-1)*v_l*eig_plus[1] + v_l*eig_plus[3] );
    flux_plus[3] = coefficient_l * ( (h_l - u_l*a_l)*eig_plus[0]
        + (gamma-1)*(u_l*u_l + v_l*v_l)*eig_plus[1] + (h_l + u_l*a_l)*eig_plus[3] );
    // Calculate right fluxes
    flux_minus[0] = coefficient_r * ( eig_minus[0]
        + 2*(gamma-1)*eig_minus[1] + eig_minus[3] );
    flux_minus[1] = coefficient_r * ( (u_r-a_r)*eig_minus[0]
        + 2*(gamma-1)*u_r*eig_minus[1] + (u_r+a_r)*eig_minus[3] );
    flux_minus[2] = coefficient_r * ( v_r*eig_minus[0]
        + 2*(gamma-1)*v_r*eig_minus[1] + v_r*eig_minus[3] );
    flux_minus[3] = coefficient_r * ( (h_r - u_r*a_r)*eig_minus[0]
        + (gamma-1)*(u_r*u_r + v_r*v_r)*eig_minus[1] + (h_r + u_r*a_r)*eig_minus[3] );

    // Add left and right fluxes to get resultant flux through face
    for (unsigned int i = 0; i < total_flux.size(); i++) {
        total_flux[i] = flux_plus[i] + flux_minus[i];
    }

    // Return result
    return total_flux;

}
