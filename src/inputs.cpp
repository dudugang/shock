#include <inputs.h>


// Constructor. This stores user inputs to the code.
Inputs::Inputs() {

    // Parameters
    dt = .000001;
    n_iterations = 100;
    n_equations = 4;
    output_rate = 1;
    mesh_file = "wedge_test.cgns";

    // Fluid properties
    gamma = 1.4;

    // Boundary conditions
    bc["in"] = BC::inflow;
    bc["out"] = BC::outflow;
    bc["sym"] = BC::wall;
    bc["wedge"] = BC::wall;

    // Inflow conditions
    rho = 1.225;
    u = 0;
    v = 0;
    p = 1e5;

    // Volume conditions
    rho_v = 1.225;
    u_v = 0;
    v_v = 0;
    p_v = 1e5;


}
