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
        double length;
        double gamma;
        int n_equations;
        int output_rate;
        string mesh_file;
        unordered_map<string, BC> bc;
};
