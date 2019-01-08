#include <vector>
#include "Algebra.h"
using std::vector;

// // // // //
// This function takes a matrix and returns a similar matrix, but with only the
// positive components. Negative components are set equal to zero.
vector<vector<double> > Algebra::return_positive(vector<vector<double> > input) {
    // Initialize
    vector<vector<double> > output = input;

    // Replace negatives with zero
    for (int i = 0; i < output.size(); i++) {
        for (int j = 0; j < output.size(); j++) {
            if (output[i][j] < 0) {
                output[i][j] = 0;
            }
        }
    }

    return output;
}

// // // // //
// This function takes a matrix and returns a similar matrix, but with only the
// negative components. Positive components are set equal to zero.
vector<vector<double> > Algebra::return_negative(vector<vector<double> > input) {
    // Initialize
    vector<vector<double> > output = input;

    // Replace positives with zero
    for (int i = 0; i < output.size(); i++) {
        for (int j = 0; j < output.size(); j++) {
            if (output[i][j] > 0) {
                output[i][j] = 0;
            }
        }
    }

    return output;
}
