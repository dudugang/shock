#include <Eigen/Dense>
#include "Algebra.h"
#include "Flux.h"
using namespace Eigen;

// // // // //
// This function calculates the flux f_i+1/2 (the flux through the right cell
// face) using Steger and Warming's Flux Vector Splitting Method. This function
// assumes a first-order reconstruction, where q is constant in a cell.
Vector3d Flux::calculate_f_right(Vector3d q, Vector3d q_right, double gamma) {
    // Calculate needed infomation from cell i
    Matrix3d r_eigvecs = calculate_right_eigenvectors(q, gamma);
    Matrix3d l_eigvecs = calculate_left_eigenvectors(q, gamma);
    Matrix3d eigenvalues = calculate_eigenvalues(q, gamma);
    Matrix3d plus_eigenvalues = Algebra::return_positive(eigenvalues);

    // Calculate needed information from cell i+1
    Matrix3d r_eigvecs_right = calculate_right_eigenvectors(q_right, gamma);
    Matrix3d l_eigvecs_right = calculate_left_eigenvectors(q_right, gamma);
    Matrix3d eigenvalues_right = calculate_eigenvalues(q_right, gamma);
    Matrix3d minus_eigenvalues_right = Algebra::return_negative(eigenvalues_right);

    // Calculate flux
    Vector3d f_right = r_eigvecs * plus_eigenvalues * l_eigvecs * q
                     + r_eigvecs_right * minus_eigenvalues_right * l_eigvecs_right * q_right;

    return f_right;
}

// // // // //
// This function calculates the flux f_i-1/2 (the flux through the left cell
// face) using Steger and Warming's Flux Vector Splitting Method. This function
// assumes a first-order reconstruction, where q is constant in a cell.
Vector3d Flux::calculate_f_left(Vector3d q, Vector3d q_left, double gamma) {
    // Calculate needed infomation from cell i
    Matrix3d r_eigvecs = calculate_right_eigenvectors(q, gamma);
    Matrix3d l_eigvecs = calculate_left_eigenvectors(q, gamma);
    Matrix3d eigenvalues = calculate_eigenvalues(q, gamma);
    Matrix3d minus_eigenvalues = Algebra::return_negative(eigenvalues);

    // Calculate needed information from cell i-1
    Matrix3d r_eigvecs_left = calculate_right_eigenvectors(q_left, gamma);
    Matrix3d l_eigvecs_left= calculate_left_eigenvectors(q_left, gamma);
    Matrix3d eigenvalues_left= calculate_eigenvalues(q_left, gamma);
    Matrix3d plus_eigenvalues_left= Algebra::return_positive(eigenvalues_left);

    // Calculate flux
    Vector3d f_left = r_eigvecs * minus_eigenvalues * l_eigvecs * q
                    + r_eigvecs_left * plus_eigenvalues_left * l_eigvecs_left * q_left;

    return f_left;
}

// // // // //
// This function calculates the matrix r_eigvecs which is the matrix of right
// eigenvectors in the columns, for the flux jacobian of the 1D Euler equations.
// The input is a q vector, which contains three components, rho, rho*u, and E.
// These equations were generated with a Matlab code, and the results were
// automatically output into C++ code format, so this was not written by hand.
Matrix3d Flux::calculate_right_eigenvectors(Vector3d q, double gamma) {
    // Initialize
    Matrix3d r_eigvecs;

    // Calculate
    r_eigvecs(0,0) = (q(0)*q(0))*1.0/(q(1)*q(1))*2.0;
    r_eigvecs(0,1) = ((q(0)*q(0))*(q(1)*q(1))*6.0-(q(0)*q(0))*(q(1)*q(1))*gamma*2.0+(q(0)*q(0)*q(0))*q(2)*gamma*4.0)/((q(1)*q(1)*q(1)*q(1))*gamma*-4.0+q(1)*q(1)*q(1)*q(1)+(q(1)*q(1)*q(1)*q(1))*(gamma*gamma)*3.0+(q(0)*q(0))*(q(2)*q(2))*(gamma*gamma)*4.0+q(0)*(q(1)*q(1))*q(2)*gamma*8.0-q(0)*(q(1)*q(1))*q(2)*(gamma*gamma)*8.0)-((q(0)*q(0))*q(1)*(q(1)*2.0+sqrt(2.0)*sqrt(gamma)*sqrt((q(0)*q(2)*2.0-q(1)*q(1))*(gamma-1.0)))*2.0)/((q(1)*q(1)*q(1)*q(1))*gamma*-4.0+q(1)*q(1)*q(1)*q(1)+(q(1)*q(1)*q(1)*q(1))*(gamma*gamma)*3.0+(q(0)*q(0))*(q(2)*q(2))*(gamma*gamma)*4.0+q(0)*(q(1)*q(1))*q(2)*gamma*8.0-q(0)*(q(1)*q(1))*q(2)*(gamma*gamma)*8.0);
    r_eigvecs(0,2) = ((q(0)*q(0))*(q(1)*q(1))*6.0-(q(0)*q(0))*(q(1)*q(1))*gamma*2.0+(q(0)*q(0)*q(0))*q(2)*gamma*4.0)/((q(1)*q(1)*q(1)*q(1))*gamma*-4.0+q(1)*q(1)*q(1)*q(1)+(q(1)*q(1)*q(1)*q(1))*(gamma*gamma)*3.0+(q(0)*q(0))*(q(2)*q(2))*(gamma*gamma)*4.0+q(0)*(q(1)*q(1))*q(2)*gamma*8.0-q(0)*(q(1)*q(1))*q(2)*(gamma*gamma)*8.0)-((q(0)*q(0))*q(1)*(q(1)*2.0-sqrt(2.0)*sqrt(gamma)*sqrt((q(0)*q(2)*2.0-q(1)*q(1))*(gamma-1.0)))*2.0)/((q(1)*q(1)*q(1)*q(1))*gamma*-4.0+q(1)*q(1)*q(1)*q(1)+(q(1)*q(1)*q(1)*q(1))*(gamma*gamma)*3.0+(q(0)*q(0))*(q(2)*q(2))*(gamma*gamma)*4.0+q(0)*(q(1)*q(1))*q(2)*gamma*8.0-q(0)*(q(1)*q(1))*q(2)*(gamma*gamma)*8.0);
    r_eigvecs(1,0) = (q(0)*2.0)/q(1);
    r_eigvecs(1,1) = -((q(1)*2.0+sqrt(2.0)*sqrt(gamma)*sqrt((q(0)*q(2)*2.0-q(1)*q(1))*(gamma-1.0)))*((q(0)*q(0))*(q(1)*q(1))+(q(0)*q(0))*(q(1)*q(1))*gamma-(q(0)*q(0)*q(0))*q(2)*gamma*2.0))/(q(0)*((q(1)*q(1)*q(1)*q(1))*gamma*-4.0+q(1)*q(1)*q(1)*q(1)+(q(1)*q(1)*q(1)*q(1))*(gamma*gamma)*3.0+(q(0)*q(0))*(q(2)*q(2))*(gamma*gamma)*4.0+q(0)*(q(1)*q(1))*q(2)*gamma*8.0-q(0)*(q(1)*q(1))*q(2)*(gamma*gamma)*8.0))+(q(0)*q(1)*(-(q(1)*q(1))*gamma+(q(1)*q(1))*2.0+(q(1)*q(1))*(gamma*gamma)+q(0)*q(2)*gamma*2.0-q(0)*q(2)*(gamma*gamma)*2.0)*2.0)/((q(1)*q(1)*q(1)*q(1))*gamma*-4.0+q(1)*q(1)*q(1)*q(1)+(q(1)*q(1)*q(1)*q(1))*(gamma*gamma)*3.0+(q(0)*q(0))*(q(2)*q(2))*(gamma*gamma)*4.0+q(0)*(q(1)*q(1))*q(2)*gamma*8.0-q(0)*(q(1)*q(1))*q(2)*(gamma*gamma)*8.0);
    r_eigvecs(1,2) = -((q(1)*2.0-sqrt(2.0)*sqrt(gamma)*sqrt((q(0)*q(2)*2.0-q(1)*q(1))*(gamma-1.0)))*((q(0)*q(0))*(q(1)*q(1))+(q(0)*q(0))*(q(1)*q(1))*gamma-(q(0)*q(0)*q(0))*q(2)*gamma*2.0))/(q(0)*((q(1)*q(1)*q(1)*q(1))*gamma*-4.0+q(1)*q(1)*q(1)*q(1)+(q(1)*q(1)*q(1)*q(1))*(gamma*gamma)*3.0+(q(0)*q(0))*(q(2)*q(2))*(gamma*gamma)*4.0+q(0)*(q(1)*q(1))*q(2)*gamma*8.0-q(0)*(q(1)*q(1))*q(2)*(gamma*gamma)*8.0))+(q(0)*q(1)*(-(q(1)*q(1))*gamma+(q(1)*q(1))*2.0+(q(1)*q(1))*(gamma*gamma)+q(0)*q(2)*gamma*2.0-q(0)*q(2)*(gamma*gamma)*2.0)*2.0)/((q(1)*q(1)*q(1)*q(1))*gamma*-4.0+q(1)*q(1)*q(1)*q(1)+(q(1)*q(1)*q(1)*q(1))*(gamma*gamma)*3.0+(q(0)*q(0))*(q(2)*q(2))*(gamma*gamma)*4.0+q(0)*(q(1)*q(1))*q(2)*gamma*8.0-q(0)*(q(1)*q(1))*q(2)*(gamma*gamma)*8.0);
    r_eigvecs(2,0) = 1.0;
    r_eigvecs(2,1) = 1.0;
    r_eigvecs(2,2) = 1.0;

    return r_eigvecs;
}

// // // // //
// This function calculates the matrix l_eigvecs which is the matrix of left
// eigenvectors in the rows, for the flux jacobian of the 1D Euler equations.
// The input is a q vector, which contains three components, rho, rho*u, and E.
// These equations were generated with a Matlab code, and the results were
// automatically output into C++ code format, so this was not written by hand.
Matrix3d Flux::calculate_left_eigenvectors(Vector3d q, double gamma) {
    // Initialize
    Matrix3d l_eigvecs;

    // Calculate
    l_eigvecs(0,0) = 1.0/(q(0)*q(0))*((q(1)*q(1))*gamma+q(1)*q(1)-q(0)*q(2)*gamma*2.0)*(1.0/2.0);
    l_eigvecs(0,1) = -q(1)/q(0);
    l_eigvecs(0,2) = 1.0;
    l_eigvecs(1,0) = (1.0/(q(0)*q(0))*((q(1)*q(1))*gamma+q(1)*q(1))*(1.0/2.0))/(gamma-1.0)-(1.0/(q(0)*q(0))*q(1)*(q(1)*2.0+sqrt(2.0)*sqrt(gamma)*sqrt((q(0)*q(2)*2.0-q(1)*q(1))*(gamma-1.0)))*(1.0/2.0))/(gamma-1.0);
    l_eigvecs(1,1) = (q(1)+sqrt(2.0)*sqrt(gamma)*sqrt((q(0)*q(2)*2.0-q(1)*q(1))*(gamma-1.0))*(1.0/2.0))/(q(0)*(gamma-1.0))-(q(1)*gamma)/(q(0)*(gamma-1.0));
    l_eigvecs(1,2) = 1.0;
    l_eigvecs(2,0) = (1.0/(q(0)*q(0))*((q(1)*q(1))*gamma+q(1)*q(1))*(1.0/2.0))/(gamma-1.0)-(1.0/(q(0)*q(0))*q(1)*(q(1)*2.0-sqrt(2.0)*sqrt(gamma)*sqrt((q(0)*q(2)*2.0-q(1)*q(1))*(gamma-1.0)))*(1.0/2.0))/(gamma-1.0);
    l_eigvecs(2,1) = (q(1)-sqrt(2.0)*sqrt(gamma)*sqrt((q(0)*q(2)*2.0-q(1)*q(1))*(gamma-1.0))*(1.0/2.0))/(q(0)*(gamma-1.0))-(q(1)*gamma)/(q(0)*(gamma-1.0));
    l_eigvecs(2,2) = 1.0;

    return l_eigvecs;
}

// // // // //
// This function calculates the matrix eigenvalues which is the matrix of
// eigenvalues on the diagonal, for the flux jacobian of the 1D Euler equations.
// The input is a q vector, which contains three components, rho, rho*u, and E.
// These equations were generated with a Matlab code, and the results were
// automatically output into C++ code format, so this was not written by hand.
Matrix3d Flux::calculate_eigenvalues(Vector3d q, double gamma) {
    // Initialize
    Matrix3d eigenvalues;

    // Calculate
    eigenvalues(0,0) = 1.0/(q(0)*q(0))*((q(1)*q(1))*gamma+q(1)*q(1)-q(0)*q(2)*gamma*2.0)*(1.0/2.0);
    eigenvalues(0,1) = -q(1)/q(0);
    eigenvalues(0,2) = 1.0;
    eigenvalues(1,0) = (1.0/(q(0)*q(0))*((q(1)*q(1))*gamma+q(1)*q(1))*(1.0/2.0))/(gamma-1.0)-(1.0/(q(0)*q(0))*q(1)*(q(1)*2.0+sqrt(2.0)*sqrt(gamma)*sqrt((q(0)*q(2)*2.0-q(1)*q(1))*(gamma-1.0)))*(1.0/2.0))/(gamma-1.0);
    eigenvalues(1,1) = (q(1)+sqrt(2.0)*sqrt(gamma)*sqrt((q(0)*q(2)*2.0-q(1)*q(1))*(gamma-1.0))*(1.0/2.0))/(q(0)*(gamma-1.0))-(q(1)*gamma)/(q(0)*(gamma-1.0));
    eigenvalues(1,2) = 1.0;
    eigenvalues(2,0) = (1.0/(q(0)*q(0))*((q(1)*q(1))*gamma+q(1)*q(1))*(1.0/2.0))/(gamma-1.0)-(1.0/(q(0)*q(0))*q(1)*(q(1)*2.0-sqrt(2.0)*sqrt(gamma)*sqrt((q(0)*q(2)*2.0-q(1)*q(1))*(gamma-1.0)))*(1.0/2.0))/(gamma-1.0);
    eigenvalues(2,1) = (q(1)-sqrt(2.0)*sqrt(gamma)*sqrt((q(0)*q(2)*2.0-q(1)*q(1))*(gamma-1.0))*(1.0/2.0))/(q(0)*(gamma-1.0))-(q(1)*gamma)/(q(0)*(gamma-1.0));
    eigenvalues(2,2) = 1.0;

    return eigenvalues;
}
