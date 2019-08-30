#pragma once
#include <time_integrator.h>


// Forward declare
class Cell;

// Class for the Forward Euler time integration method
class ForwardEuler : public TimeIntegrator {

    public:
        ForwardEuler();
        virtual void integrate(Cell*, double);

};
