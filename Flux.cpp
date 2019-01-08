#include <vector>
#include "Algebra.h"
#include "Flux.h"
using std::vector;

vector<double> Flux::calculate_f_right(vector<double>& q, vector<double>& q_right, double gamma) {
    // Calculate needed infomation from cell i
    vector<vector<double> > r_eigvecs = calculate_right_eigenvectors(q, gamma);
    vector<vector<double> > l_eigvecs = calculate_left_eigenvectors(q, gamma);
    vector<vector<double> > eigenvalues = calculate_eigenvalues(q, gamma);
    vector<vector<double> > plus_eigenvalues = return_positive(eigenvalues);

    // Calculate needed information from cell i+1
    vector<vector<double> > r_eigvecs_right = calculate_right_eigenvectors(q_right, gamma);
    vector<vector<double> > l_eigvecs_right = calculate_left_eigenvectors(q_right, gamma);
    vector<vector<double> > eigenvalues_right = calculate_eigenvalues(q_right, gamma);
    vector<vector<double> > minus_eigenvalues_right = return_negative(eigenvalues_right);

    // Calculate flux
    vector<double> f_right =


// // // // //
// This function calculates the matrix r_eigvecs which is the matrix of right
// eigenvectors in the columns, for the flux jacobian of the 1D Euler equations.
// The input is a q vector, which contains three components, rho, rho*u, and E.
// These equations were generated with a Matlab code, and the results were
// automatically output into C++ code format, so this was not written by hand.
vector<vector<double> > Flux::calculate_right_eigenvectors(vector<double>& q, double gamma) {
    // Initialize
    vector<vector<double> > r_eigvecs;

    // Calculate
    r_eigvecs[0][0] = (q[0]*q[0])*1.0/(q[1]*q[1])*2.0;
    r_eigvecs[0][1] = ((q[0]*q[0])*(q[1]*q[1])*6.0-(q[0]*q[0])*(q[1]*q[1])*gamma*2.0+(q[0]*q[0]*q[0])*q[2]*gamma*4.0)/((q[1]*q[1]*q[1]*q[1])*gamma*-4.0+q[1]*q[1]*q[1]*q[1]+(q[1]*q[1]*q[1]*q[1])*(gamma*gamma)*3.0+(q[0]*q[0])*(q[2]*q[2])*(gamma*gamma)*4.0+q[0]*(q[1]*q[1])*q[2]*gamma*8.0-q[0]*(q[1]*q[1])*q[2]*(gamma*gamma)*8.0)-((q[0]*q[0])*q[1]*(q[1]*2.0+sqrt(2.0)*sqrt(gamma)*sqrt((q[0]*q[2]*2.0-q[1]*q[1])*(gamma-1.0)))*2.0)/((q[1]*q[1]*q[1]*q[1])*gamma*-4.0+q[1]*q[1]*q[1]*q[1]+(q[1]*q[1]*q[1]*q[1])*(gamma*gamma)*3.0+(q[0]*q[0])*(q[2]*q[2])*(gamma*gamma)*4.0+q[0]*(q[1]*q[1])*q[2]*gamma*8.0-q[0]*(q[1]*q[1])*q[2]*(gamma*gamma)*8.0);
    r_eigvecs[0][2] = ((q[0]*q[0])*(q[1]*q[1])*6.0-(q[0]*q[0])*(q[1]*q[1])*gamma*2.0+(q[0]*q[0]*q[0])*q[2]*gamma*4.0)/((q[1]*q[1]*q[1]*q[1])*gamma*-4.0+q[1]*q[1]*q[1]*q[1]+(q[1]*q[1]*q[1]*q[1])*(gamma*gamma)*3.0+(q[0]*q[0])*(q[2]*q[2])*(gamma*gamma)*4.0+q[0]*(q[1]*q[1])*q[2]*gamma*8.0-q[0]*(q[1]*q[1])*q[2]*(gamma*gamma)*8.0)-((q[0]*q[0])*q[1]*(q[1]*2.0-sqrt(2.0)*sqrt(gamma)*sqrt((q[0]*q[2]*2.0-q[1]*q[1])*(gamma-1.0)))*2.0)/((q[1]*q[1]*q[1]*q[1])*gamma*-4.0+q[1]*q[1]*q[1]*q[1]+(q[1]*q[1]*q[1]*q[1])*(gamma*gamma)*3.0+(q[0]*q[0])*(q[2]*q[2])*(gamma*gamma)*4.0+q[0]*(q[1]*q[1])*q[2]*gamma*8.0-q[0]*(q[1]*q[1])*q[2]*(gamma*gamma)*8.0);
    r_eigvecs[1][0] = (q[0]*2.0)/q[1];
    r_eigvecs[1][1] = -((q[1]*2.0+sqrt(2.0)*sqrt(gamma)*sqrt((q[0]*q[2]*2.0-q[1]*q[1])*(gamma-1.0)))*((q[0]*q[0])*(q[1]*q[1])+(q[0]*q[0])*(q[1]*q[1])*gamma-(q[0]*q[0]*q[0])*q[2]*gamma*2.0))/(q[0]*((q[1]*q[1]*q[1]*q[1])*gamma*-4.0+q[1]*q[1]*q[1]*q[1]+(q[1]*q[1]*q[1]*q[1])*(gamma*gamma)*3.0+(q[0]*q[0])*(q[2]*q[2])*(gamma*gamma)*4.0+q[0]*(q[1]*q[1])*q[2]*gamma*8.0-q[0]*(q[1]*q[1])*q[2]*(gamma*gamma)*8.0))+(q[0]*q[1]*(-(q[1]*q[1])*gamma+(q[1]*q[1])*2.0+(q[1]*q[1])*(gamma*gamma)+q[0]*q[2]*gamma*2.0-q[0]*q[2]*(gamma*gamma)*2.0)*2.0)/((q[1]*q[1]*q[1]*q[1])*gamma*-4.0+q[1]*q[1]*q[1]*q[1]+(q[1]*q[1]*q[1]*q[1])*(gamma*gamma)*3.0+(q[0]*q[0])*(q[2]*q[2])*(gamma*gamma)*4.0+q[0]*(q[1]*q[1])*q[2]*gamma*8.0-q[0]*(q[1]*q[1])*q[2]*(gamma*gamma)*8.0);
    r_eigvecs[1][2] = -((q[1]*2.0-sqrt(2.0)*sqrt(gamma)*sqrt((q[0]*q[2]*2.0-q[1]*q[1])*(gamma-1.0)))*((q[0]*q[0])*(q[1]*q[1])+(q[0]*q[0])*(q[1]*q[1])*gamma-(q[0]*q[0]*q[0])*q[2]*gamma*2.0))/(q[0]*((q[1]*q[1]*q[1]*q[1])*gamma*-4.0+q[1]*q[1]*q[1]*q[1]+(q[1]*q[1]*q[1]*q[1])*(gamma*gamma)*3.0+(q[0]*q[0])*(q[2]*q[2])*(gamma*gamma)*4.0+q[0]*(q[1]*q[1])*q[2]*gamma*8.0-q[0]*(q[1]*q[1])*q[2]*(gamma*gamma)*8.0))+(q[0]*q[1]*(-(q[1]*q[1])*gamma+(q[1]*q[1])*2.0+(q[1]*q[1])*(gamma*gamma)+q[0]*q[2]*gamma*2.0-q[0]*q[2]*(gamma*gamma)*2.0)*2.0)/((q[1]*q[1]*q[1]*q[1])*gamma*-4.0+q[1]*q[1]*q[1]*q[1]+(q[1]*q[1]*q[1]*q[1])*(gamma*gamma)*3.0+(q[0]*q[0])*(q[2]*q[2])*(gamma*gamma)*4.0+q[0]*(q[1]*q[1])*q[2]*gamma*8.0-q[0]*(q[1]*q[1])*q[2]*(gamma*gamma)*8.0);
    r_eigvecs[2][0] = 1.0;
    r_eigvecs[2][1] = 1.0;
    r_eigvecs[2][2] = 1.0;

    return r_eigvecs;
}

// // // // //
// This function calculates the matrix l_eigvecs which is the matrix of left
// eigenvectors in the rows, for the flux jacobian of the 1D Euler equations.
// The input is a q vector, which contains three components, rho, rho*u, and E.
// These equations were generated with a Matlab code, and the results were
// automatically output into C++ code format, so this was not written by hand.
vector<vector<double> > Flux::calculate_left_eigenvectors(vector<double>& q, double gamma) {
    // Initialize
    vector<vector<double> > l_eigvecs;

    // Calculate
    l_eigvecs[0][0] = 1.0/(q[0]*q[0])*((q[1]*q[1])*gamma+q[1]*q[1]-q[0]*q[2]*gamma*2.0)*(1.0/2.0);
    l_eigvecs[0][1] = -q[1]/q[0];
    l_eigvecs[0][2] = 1.0;
    l_eigvecs[1][0] = (1.0/(q[0]*q[0])*((q[1]*q[1])*gamma+q[1]*q[1])*(1.0/2.0))/(gamma-1.0)-(1.0/(q[0]*q[0])*q[1]*(q[1]*2.0+sqrt(2.0)*sqrt(gamma)*sqrt((q[0]*q[2]*2.0-q[1]*q[1])*(gamma-1.0)))*(1.0/2.0))/(gamma-1.0);
    l_eigvecs[1][1] = (q[1]+sqrt(2.0)*sqrt(gamma)*sqrt((q[0]*q[2]*2.0-q[1]*q[1])*(gamma-1.0))*(1.0/2.0))/(q[0]*(gamma-1.0))-(q[1]*gamma)/(q[0]*(gamma-1.0));
    l_eigvecs[1][2] = 1.0;
    l_eigvecs[2][0] = (1.0/(q[0]*q[0])*((q[1]*q[1])*gamma+q[1]*q[1])*(1.0/2.0))/(gamma-1.0)-(1.0/(q[0]*q[0])*q[1]*(q[1]*2.0-sqrt(2.0)*sqrt(gamma)*sqrt((q[0]*q[2]*2.0-q[1]*q[1])*(gamma-1.0)))*(1.0/2.0))/(gamma-1.0);
    l_eigvecs[2][1] = (q[1]-sqrt(2.0)*sqrt(gamma)*sqrt((q[0]*q[2]*2.0-q[1]*q[1])*(gamma-1.0))*(1.0/2.0))/(q[0]*(gamma-1.0))-(q[1]*gamma)/(q[0]*(gamma-1.0));
    l_eigvecs[2][2] = 1.0;

    return l_eigvecs;
}

// // // // //
// This function calculates the matrix eigenvalues which is the matrix of
// eigenvalues on the diagonal, for the flux jacobian of the 1D Euler equations.
// The input is a q vector, which contains three components, rho, rho*u, and E.
// These equations were generated with a Matlab code, and the results were
// automatically output into C++ code format, so this was not written by hand.
vector<vector<double> > Flux::calculate_eigenvalues(vector<double>& q, double gamma) {
    // Initialize
    vector<vector<double> > eigenvalues;
    eigenvalues.resize(3);
    eigenvalues[0].resize(3);
    eigenvalues[1].resize(3);
    eigenvalues[2].resize(3);

    // Calculate
    eigenvalues[0][0] = 1.0/(q[0]*q[0])*((q[1]*q[1])*gamma+q[1]*q[1]-q[0]*q[2]*gamma*2.0)*(1.0/2.0);
    eigenvalues[0][1] = -q[1]/q[0];
    eigenvalues[0][2] = 1.0;
    eigenvalues[1][0] = (1.0/(q[0]*q[0])*((q[1]*q[1])*gamma+q[1]*q[1])*(1.0/2.0))/(gamma-1.0)-(1.0/(q[0]*q[0])*q[1]*(q[1]*2.0+sqrt(2.0)*sqrt(gamma)*sqrt((q[0]*q[2]*2.0-q[1]*q[1])*(gamma-1.0)))*(1.0/2.0))/(gamma-1.0);
    eigenvalues[1][1] = (q[1]+sqrt(2.0)*sqrt(gamma)*sqrt((q[0]*q[2]*2.0-q[1]*q[1])*(gamma-1.0))*(1.0/2.0))/(q[0]*(gamma-1.0))-(q[1]*gamma)/(q[0]*(gamma-1.0));
    eigenvalues[1][2] = 1.0;
    eigenvalues[2][0] = (1.0/(q[0]*q[0])*((q[1]*q[1])*gamma+q[1]*q[1])*(1.0/2.0))/(gamma-1.0)-(1.0/(q[0]*q[0])*q[1]*(q[1]*2.0-sqrt(2.0)*sqrt(gamma)*sqrt((q[0]*q[2]*2.0-q[1]*q[1])*(gamma-1.0)))*(1.0/2.0))/(gamma-1.0);
    eigenvalues[2][1] = (q[1]-sqrt(2.0)*sqrt(gamma)*sqrt((q[0]*q[2]*2.0-q[1]*q[1])*(gamma-1.0))*(1.0/2.0))/(q[0]*(gamma-1.0))-(q[1]*gamma)/(q[0]*(gamma-1.0));
    eigenvalues[2][2] = 1.0;

    return eigenvalues;
}
