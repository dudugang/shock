#include <iostream>
#include <output.h>
using std::cout;
using std::endl;

int main(int argc, char* argv[]) {
    // Get inputs
    Inputs inputs;

    // Initialize flowfield
    Flowfield flow(inputs);

    // Main loop: run flux function and time integrator for every iteration
    for (int i = 0; i < inputs.n_iterations; i++) {

        // Calculate fluxes
        flow.calculate_flux();

        // Integrate in time
        flow.apply_time_integrator();

    }

    // Output results
    Output output;
    output.print(flow);

    return 0;
}
