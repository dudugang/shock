#pragma once
#include <vector>
using std::vector;


// Forward declare
class Cell;

// Class for integrating the semi-discrete form of the governing equations in
// time.
class TimeIntegrator {

    public:
        static void initialize();
        static void forward_euler(Cell*, double);

    private:
        static vector<double> flux_integral;

};
