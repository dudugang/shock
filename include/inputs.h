#include <iostream>
#include <string>
#include <vector>
#include <H5Cpp.h>
using namespace H5;
using std::cout;
using std::endl;
using std::string;
using std::vector;

// Class for reading/storing user inputs to the code.
class Inputs {
    public:
        Inputs();
        void read_mesh();
        double dt;
        int n_cells;
        int n_iterations;
        double length;
        double dx;
        double gamma;
        int n_equations;
        int output_rate;
        vector<double> q_left;
        vector<double> q_right;
};
