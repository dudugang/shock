#include <inputs.h>


// Constructor. This stores user inputs to the code.
Inputs::Inputs() {

    // Parameters
    dt = .00001;
    n_iterations = 2000;
    n_equations = 4;
    output_rate = 20;
    case_file = "box.cgns";

    // Fluid properties
    gamma = 1.4;

    // Debug
    print_id = 1;

    // Boundary conditions
    bc["inlet"]   = BC::inflow;
    bc["outflow"] = BC::outflow;
    bc["in"]   = BC::inflow;
    bc["out"]  = BC::outflow;
    bc["sym"]  = BC::wall;
    bc["wall"] = BC::wall;

    // Volume conditions
    vc["interior"] = {1, 0, 0, 1e5};
    vc["driver"] = {1, 0, 0, 1e5};
    vc["driven"] = {1, 0, 0, 1e5};

    // Inflow conditions
    rho = 1;
    u = 1200;
    v = 0;
    p = 1e5;

}
