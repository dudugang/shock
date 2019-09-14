#include <face.h>
#include <cell.h>
#include <ghost.h>


// Constructor
Face::Face(Point point1, Point point2) {

    // Make sure that the points are oriented according to the sign conventions
    // illustrated in the header file
    if (point1.x > point2.x or point1.y > point2.y) {
        this->point1 = point1;
        this->point2 = point2;
    } else {
        this->point1 = point2;
        this->point2 = point1;
    }

    flux.resize(4);
    q_left.resize(4);
    q_right.resize(4);
    neighbors.reserve(2);

    // Find centroid of face
    center = Geometry::find_midpoint(point1, point2);

}


// Calculate the angle of the normal vector for a face, making sure that the
// vector points towards the cell with the higher cell ID.
void Face::find_normal_vector(unordered_map<int, Cell*> &cells,
    unordered_map<int, Ghost*> &ghosts, int n_cells) {

    // Get smaller of the two neighbor cell ID's
    int id = neighbors[0];

    // Get cell center of whichever cell has the lowest cell ID. The reason why
    // lowest was chosen instead of highest is because ghost cells always have a
    // higher cell ID than flow cells, which would be a problem since the
    // geometry of a ghost cell is never defined, since it doesn't need to be.
    Point cell_center;
    if (id <= n_cells) {
        cell_center = cells[id]->center;
    } else {
        cell_center = ghosts[id]->center;
    }

    // Find area. This is just the length of the line segment between the two
    // points, since this code is 2D (for now).
    double dx = point1.x - point2.x;
    double dy = point1.y - point2.y;
    area = std::sqrt(dx*dx + dy*dy);

    // Find one of the possible normal vector angles
    double theta1 = std::atan2(dy, dx) + Geometry::pi/2;

    // Find other possible normal vector angle
    double theta2 = theta1 + Geometry::pi;
    // Make sure this doesn't equal or exceed 2pi
    if (theta2 >= 2*Geometry::pi) {
        theta2 -= 2*Geometry::pi;
    }

    // Get angle between x-axis and the line connecting the chosen cell center
    // and the face center
    double theta_cell = std::atan2(center.y - cell_center.y, center.x - cell_center.x);
    // Make sure it's not negative
    if (theta_cell < 0) {
        theta_cell += 2*Geometry::pi;
    }

    // Choose the theta that gives a normal vector pointing towards the cell
    // with larger cell ID. This is done by choosing the theta which minimizes
    // the difference between the normal vector angle and the angle between the
    // cell and face with respect to the x-axis.
    if (Geometry::angle_between(theta_cell, theta1)
            < Geometry::angle_between(theta_cell, theta2)) {
        theta = theta1;
    } else {
        theta = theta2;
    }

    // Store sin(theta) and cos(theta) for convenience
    sintheta = std::sin(theta);
    costheta = std::cos(theta);

}


// Sort neighbor ID's from small to large
void Face::sort_neighbors() { std::sort(neighbors.begin(), neighbors.end()); }


// Check if face nodes contain certain points
bool Face::contains(Point p1, Point p2) {

    // Output
    bool result;

    // Check if object contains the points in any order
    if ((point1 == p1) and (point2 == p2)) {
        result = true;
    } else if ((point1 == p2) and (point2 == p1)) {
        result = true;
    } else {
        result = false;
    }

    return result;

}


// Check if two faces equal by comparing midpoints
bool operator== (const Face &face1, const Face &face2) {
    return (face1.center.x == face2.center.x
            and face1.center.y == face2.center.y);
}
