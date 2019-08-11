#include <iostream>
#include <mesh_reader.h>
#include <output.h>
using std::cout;
using std::endl;

// This is the main, the highest level of the code. This function glues together
// the rest of the code, without getting into implementation details.
int main(int argc, char* argv[]) {

    // Get inputs
    Inputs inputs;

    // Read in mesh
    MeshReader mesh_reader(inputs.mesh_file);
    mesh_reader.create_mesh();

    // Initialize flowfield
    Flowfield flow(inputs, mesh_reader);

    // Initialize output
    Output output(inputs, flow.n_cells);

    // Initialize fluxes
    Flux flux;

    // Main loop: run flux function and time integrator for every iteration
    for (int i = 0; i < inputs.n_iterations; i++) {

        // Calculate fluxes
        flow.calculate_flux(flux);

        // Reconstruct flowfield
        flow.apply_reconstruction();

        // Integrate in time
        flow.apply_time_integrator();
/*

        // Write to file
        output.write(flow, i);

        // Output information to stdout
        output.print(flow, i);

*/
    }

    // Output results
    //output.final_print(flow);
    return 0;
}
