#include <Eigen/Dense>
using namespace Eigen;

class Flux {
    public:
        Matrix3d calculate_right_eigenvectors(Vector3d, double);
        Matrix3d calculate_left_eigenvectors (Vector3d, double);
        Matrix3d calculate_eigenvalues       (Vector3d, double);
        Vector3d calculate_f_right(Vector3d, Vector3d, double);
        Vector3d calculate_f_left (Vector3d, Vector3d, double);
};
