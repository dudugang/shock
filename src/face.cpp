#include <face.h>
using std::cout;
using std::endl;

Face::Face(int face_id, Volume *left_volume, Volume *right_volume) {
    cout << "Creating face with ID " << face_id << endl;
    this->face_id = face_id;
    this->left_volume = left_volume;
    this->right_volume = right_volume;
    flux.resize(3);
}
