#include <volume.h>


// Constructor
Volume::Volume(vector<Point> vertices, vector<double> q, vector<Face*> faces,
    int volume_id) {
    this->vertices = vertices;
    this->q = q;
    this->faces = faces;
    this->volume_id = volume_id;
}
