#include <iostream>
#include <flowfield.h>
#include <flux.h>
#include <inputs.h>
#include <mesh_reader.h>
#include <output.h>
#include <time_integrator.h>
using std::cout;
using std::endl;

// This is the main, the highest level of the code. This function glues together
// the rest of the code, without getting into implementation details.
int main() {

    // Get inputs
    Inputs inputs;

    // Read in mesh
    MeshReader mesh_reader(inputs.case_file);
    mesh_reader.create_mesh(inputs);

    // Initialize flowfield
    Flowfield flow(inputs, mesh_reader);

    // Initialize fluxes
    Flux flux;

    // Initialize time integrator
    TimeIntegrator* time_integrator = TimeIntegrator::choose_time_integrator(inputs);

    // Initialize output
    Output output(inputs, flow);

    // Main loop: run flux function and time integrator for every iteration
    for (int i = 1; i <= inputs.n_iterations; i++) {

        // Reconstruct flowfield
        flow.apply_reconstruction();

        // Calculate fluxes
        flow.calculate_flux(flux);

        // Integrate in time
        flow.apply_time_integrator(time_integrator);

        // Write to file
        output.write(flow, i);

        // Output information to stdout and save time
        output.print(flow, i);
        output.add_time(flow.time);

    }

    // Output results
    //output.write_results(inputs.case_file, flow.cells, flow.n_cells);
    return 0;
}
