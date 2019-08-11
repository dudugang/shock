#include <cell.h>


// Constructor
Cell::Cell(vector<Point> vertices, vector<double> q, vector<Face*> faces,
    int volume_id) : Volume(vertices, q, faces, volume_id) {
    this->type = "flow";
    cout << "Creating cell with ID " << volume_id << endl;

    // Calculate cell center, assuming quadrilateral cells
    // TODO: Stop assuming quadrilateral cells
    center = Geometry::find_quad_centroid(vertices[0], vertices[1], vertices[2],
        vertices[3]);

    // Calculate cell volume
    volume = Geometry::find_quad_area(vertices[0], vertices[1], vertices[2],
        vertices[3]);

}

void Cell::update() {}
