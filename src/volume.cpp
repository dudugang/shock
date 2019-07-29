#include <volume.h>
#include <vector>
using std::cout;
using std::endl;
using std::vector;

// Constructor
Volume::Volume(vector<double> center, vector<double> q, int volume_id,
    vector<Volume*> neighbors, Face *left_face, Face *right_face) {
    cout << "Creating volume with ID " << volume_id << endl;
    this->center = center;
    this->q = q;
    this->volume_id = volume_id;
    this->neighbors = neighbors;
    this->left_face = left_face;
    this->right_face = right_face;
}
