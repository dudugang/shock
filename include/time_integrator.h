#pragma once
#include <vector>
#include <inputs.h>
using std::vector;


// Forward declare
class Cell;

// Class for integrating the semi-discrete form of the governing equations in
// time.
class TimeIntegrator {

    public:
        static TimeIntegrator* choose_time_integrator(Inputs);
        virtual void integrate(Cell*, double) = 0;

    protected:
        vector<double> flux_integral;

};
