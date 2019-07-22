#include <inputs.h>

// Constructor. This stores user inputs to the code.
Inputs::Inputs() {

    // Parameters
    dt = .0001;
    n_cells = 100;
    n_iterations = 5000;
    length = 1;
    dx = length/n_cells;
    gamma = 1.4;

    // Initial conditions of conserved variables rho, rho*u, and rho*e
    q_left[0]  = 1.5;
    q_left[1]  = 0.00001;
    q_left[2]  = 3.7e6;
    q_right[0] = 1;
    q_right[1] = 0.00001;
    q_right[2] = 2.5e6;
}
