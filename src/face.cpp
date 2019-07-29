#include <face.h>
using std::cout;
using std::endl;


// Constructor
Face::Face(int face_id, Volume *left_volume, Volume *right_volume,
    vector<double> point1, vector<double> point2) {
    cout << "Creating face with ID " << face_id << endl;
    this->face_id = face_id;
    this->left_volume = left_volume;
    this->right_volume = right_volume;
    flux.resize(3);
    this->point1 = point1;
    this->point2 = point2;

    // Calculate angle of normal vector, using sign conventions illustrated in
    // the header file for this class
    double dx = point1[0] - point2[0];
    double dy = point1[1] - point2[1];
    ds = std::sqrt(dx*dx + dy*dy);
    sintheta = dy/ds;
    costheta = dx/ds;
    theta = std::atan2(sintheta, costheta);
}
