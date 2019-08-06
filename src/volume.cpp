#include <volume.h>
#include <vector>
using std::cout;
using std::endl;
using std::vector;

// Constructor
Volume::Volume(vector<Point> vertices, vector<double> q, vector<Face*> faces,
    int volume_id) {
    cout << "Creating volume with ID " << volume_id << endl;
    this->vertices = vertices;
    this->q = q;
    this->faces = faces;
    this->volume_id = volume_id;
}
