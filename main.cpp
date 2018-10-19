#include "Flow.h"

int main(int argc, char* argv[]) {
    // Make instance of Flow class and initialize inputs
    Flow flow;

    // Generates grid and initial conditions from inputs
    flow.initialize();

    // Run computations
    flow.solve();

    // Output solution
    flow.output();

    return 0;
}
