#include <inputs.h>


// Constructor. This stores user inputs to the code.
Inputs::Inputs() {

    // Parameters
    dt = .0001;
    n_iterations = 1;
    n_equations = 4;
    output_rate = 1;
    case_file = "wedge_test.cgns";

    // Fluid properties
    gamma = 1.4;

    // Boundary conditions
    bc["in"] = BC::inflow;
    bc["out"] = BC::outflow;
    bc["sym"] = BC::wall;
    bc["wedge"] = BC::wall;

    // Inflow conditions
    rho = .1;
    u = 100;
    v = 0;
    p = 1e4;

    // Volume conditions
    rho_v = .1;
    u_v = 0;
    v_v = 0;
    p_v = 1e4;


}
