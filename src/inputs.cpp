#include <inputs.h>


// Constructor. This stores user inputs to the code.
Inputs::Inputs() {

    // Parameters
    dt = .000001;
    n_iterations = 1;
    length = 1;
    gamma = 1.4;
    n_equations = 4;
    output_rate = 1;
    mesh_file = "wedge_test.cgns";

    // Boundary conditions
    //bc["in"] = BC::inflow;
    //bc["out"] = BC::inflow;
    //bc["sym"] = BC::inflow;
    //bc["wedge"] = BC::inflow;
    bc["wall"] = BC::inflow;

    // Initial conditions of conserved variables rho, rho*u, and rho*e
    q_left.resize(4);
    q_right.resize(4);
    q_left[0]  = 1.5;
    q_left[1]  = 0.00001;
    q_left[2]  = 0.00001;
    q_left[3]  = 3.7e6;
    q_right[0] = 1;
    q_right[1] = 0.00001;
    q_right[2] = 0.00001;
    q_right[3] = 2.5e6;
}
