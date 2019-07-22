#include <vector>
using std::vector;

// Class for reading/storing user inputs to the code.
class Inputs {
    public:
        Inputs();
        double dt;
        int n_cells;
        int n_iterations;
        double length;
        double dx;
        double gamma;
        vector<double> q_left;
        vector<double> q_right;
};
