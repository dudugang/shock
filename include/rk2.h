#pragma once
#include <unordered_map>
#include <flowfield.h>
#include <time_integrator.h>
using std::unordered_map;


// Forward declare
class Cell;

// Class for the Forward Euler time integration method
class RK2 : public TimeIntegrator {

    public:
        RK2();
        virtual void integrate(Flowfield&, unordered_map<int, Cell*>&, Flux& flux, double);

};
