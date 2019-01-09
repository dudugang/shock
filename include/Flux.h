#include <Eigen/Dense>
using namespace Eigen;

class Flux {
    private:
        static Matrix3d calculate_right_eigenvectors(Vector3d, double);
        static Matrix3d calculate_left_eigenvectors (Vector3d, double);
        static Matrix3d calculate_eigenvalues       (Vector3d, double);
    public:
        static Vector3d calculate_f_right(Vector3d, Vector3d, double);
        static Vector3d calculate_f_left (Vector3d, Vector3d, double);
};
