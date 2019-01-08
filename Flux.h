#include <vector>
using std::vector;

class Flux {
    public:
        vector<vector<double> > calculate_right_eigenvectors(vector<double>&, double);
        vector<vector<double> > calculate_left_eigenvectors (vector<double>&, double);
        vector<vector<double> > calculate_eigenvalues       (vector<double>&, double);
        vector<double> calculate_f_right(vector<double>&, vector<double>&, double);
        vector<double> calculate_f_left(vector<double>&, vector<double>&, double);
};
