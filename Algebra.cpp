#include <Eigen/Dense>
#include "Algebra.h"
using namespace Eigen;

// // // // //
// This function takes a matrix and returns a similar matrix, but with only the
// positive components. Negative components are set equal to zero.
MatrixXd Algebra::return_positive(MatrixXd input) {
    // Initialize
    MatrixXd output = input;

    // Replace negatives with zero
    for (int i = 0; i < output.rows(); i++) {
        for (int j = 0; j < output.cols(); j++) {
            if (output(i,j) < 0) {
                output(i,j) = 0;
            }
        }
    }

    return output;
}

// // // // //
// This function takes a matrix and returns a similar matrix, but with only the
// negative components. Positive components are set equal to zero.
MatrixXd Algebra::return_negative(MatrixXd input) {
    // Initialize
    MatrixXd output = input;

    // Replace negatives with zero
    for (int i = 0; i < output.rows(); i++) {
        for (int j = 0; j < output.cols(); j++) {
            if (output(i,j) > 0) {
                output(i,j) = 0;
            }
        }
    }

    return output;
}
