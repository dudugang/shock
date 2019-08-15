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
    bc["in"] = BC::inflow;
    bc["out"] = BC::inflow;
    bc["sym"] = BC::inflow;
    bc["wedge"] = BC::inflow;
}
