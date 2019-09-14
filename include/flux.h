#pragma once
#include <cmath>
#include <iostream>
#include <vector>
using std::cout;
using std::endl;
using std::sqrt;
using std::vector;

// Class containing static methods for computing fluxes between cells.
class Flux {
    public:
        Flux();
        vector<double> steger_warming(vector<double>, vector<double>, double);
        vector<double> eig_plus;
        vector<double> eig_minus;
        vector<double> flux_plus;
        vector<double> flux_minus;
        vector<double> total_flux;
    private:
        vector<double> eigvalues_l;
        vector<double> eigvalues_r;
};
