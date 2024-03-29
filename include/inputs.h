#pragma once
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <bc.h>
using std::cout;
using std::endl;
using std::string;
using std::unordered_map;
using std::vector;

// Class for reading/storing user inputs to the code.
class Inputs {
    public:
        Inputs();
        void read_mesh();
        double dt;
        int n_iterations;
        int n_equations;
        int output_rate;
        string case_file;
        double gamma;
        string time_integrator;
        int print_id;
        unordered_map<string, BC> bc;
        unordered_map<string, vector<double> > vc;
        double rho;
        double u;
        double v;
        double p;
        double rho_v;
        double u_v;
        double v_v;
        double p_v;
};
