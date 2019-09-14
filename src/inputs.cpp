#include <inputs.h>


// Constructor. This stores user inputs to the code.
Inputs::Inputs() {

    // Parameters
    dt = 1e-5;
    n_iterations = 1000;
    n_equations = 4;
    output_rate = 20;
    case_file = "wedge_30deg.cgns";

    // Fluid properties
    gamma = 1.4;

    // Time integrator
    time_integrator = "ForwardEuler";

    // Debug
    print_id = 1;

    // Boundary conditions
    bc["inlet"]   = BC::inflow;
    bc["outflow"] = BC::outflow;
    bc["outlet"] = BC::outflow;
    bc["in"]    = BC::inflow;
    bc["out"]   = BC::outflow;
    bc["sym"]   = BC::wall;
    bc["wall"]  = BC::wall;
    bc["wedge"] = BC::wall;

    // Volume conditions
    vc["interior"] = {1, 1000, 0, 1e5};
    vc["driver"]   = {1, 1000, 0, 1e5};
    vc["driven"]   = {1, 1000, 0, 1e5};

    // Inflow conditions
    rho = 1;
    u = 1000;
    v = 0;
    p = 1e5;

}
