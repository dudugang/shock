#include <volume.h>
#include <face.h>


// Constructor
Volume::Volume(vector<Point> vertices, vector<double> q, vector<Face*> faces,
    int id) {
    this->vertices = vertices;
    this->q = q;
    this->faces = faces;
    this->id = id;
}
