#include <cell.h>
#include <face.h>


// Constructor
Cell::Cell(vector<Point> vertices, vector<double> q, vector<Face*> faces,
    int id) : Volume(vertices, q, faces, id) {
    this->type = "flow";
    cout << "Creating cell with ID " << id << endl;

    // Calculate cell center, assuming quadrilateral cells
    // TODO: Stop assuming quadrilateral cells
    center = Geometry::find_quad_centroid(vertices[0], vertices[1], vertices[2],
        vertices[3]);

    // Calculate cell volume
    volume = Geometry::find_quad_area(vertices[0], vertices[1], vertices[2],
        vertices[3]);

}
