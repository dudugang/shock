#include <cell.h>


// Constructor
Cell::Cell(vector<Point> vertices, vector<double> q, vector<Face*> faces,
    int volume_id) : Volume(vertices, q, faces, volume_id) {
    this->type = "flow";
    cout << "Creating cell with ID " << volume_id << endl;
}

void Cell::update() {}
