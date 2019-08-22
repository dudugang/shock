#include <inputs.h>


// Constructor. This stores user inputs to the code.
Inputs::Inputs() {

    // Parameters
    dt = .000005;
    n_iterations = 1000;
    n_equations = 4;
    output_rate = 10;
    case_file = "shock_tube.cgns";

    // Fluid properties
    gamma = 1.4;

    // Boundary conditions
    bc["in"]   = BC::inflow;
    bc["out"]  = BC::outflow;
    bc["sym"]  = BC::wall;
    bc["wall"] = BC::outflow;

    // Volume conditions
    vc["driver"] = {1, 0, 0, 1e5};
    vc["driven"] = {1, 0, 0, 1e4};

    // Inflow conditions
    rho = 1;
    u = 1200;
    v = 0;
    p = 1e5;

}
