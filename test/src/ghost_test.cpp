#include <cmath>
#include <unordered_map>
#include <vector>
#include <cell.h>
#include <face.h>
#include <geometry.h>
#include <ghost.h>
#include <point.h>
#include <gtest/gtest.h>
using std::unordered_map;
using std::vector;


// Function for creating test data that is used multiple times
void create_ghost(vector<double> q_neighbor, Ghost* &ghost,
    unordered_map<int, Cell*> &cells) {

    // Create face and place in vector of faces
    Point a(0, 0);
    Point b(1, 0);
    Face* face = new Face(a, b);
    vector<Face*> faces = {face};

    // Set angle of face
    double theta = Geometry::pi / 2;
    face->sintheta = std::sin(theta);
    face->costheta = std::cos(theta);

    // Create neighbor cell and place in map of cells
    Point c(1, 1);
    Point d(0, 1);
    vector<Point> vertices = {a, b, c, d};
    int neighbor_id = 0;
    Cell* cell = new Cell(vertices, q_neighbor, faces, neighbor_id);
    cells[neighbor_id] = cell;

    // Create ghost cell
    vector<double> q;
    q.resize(4);
    int id = 1;
    ghost = new Ghost(vertices, q, faces, id);

    // Add these cells to face's neighbors
    face->neighbors.push_back(neighbor_id);
    face->neighbors.push_back(id);

}


TEST(update_wall_test, ShouldDoNothingForTangentialFlow) {

    // Create ghost object
    vector<double> q_neighbor = {1, 1, 0, 1e5};
    Ghost* ghost;
    unordered_map<int, Cell*> cells;
    create_ghost(q_neighbor, ghost, cells);

    // Need gamma for this BC
    double gamma = 1.4;

    // Apply boundary condition
    ghost->update_wall(cells, gamma);
    vector<double> ghost_q = ghost->q;
    vector<double> ghost_q_correct = {1, 1, 0, 1e5};

    // Test
    double error = 1e-8;
    for (size_t i = 0; i < ghost_q.size(); i++) {
        ASSERT_NEAR(ghost_q_correct[i], ghost_q[i], error);
    }

}

TEST(update_wall_test, ShouldFlipSignForNormalFlow) {

    // Create ghost object
    vector<double> q_neighbor = {1, 0, 1, 1e5};
    Ghost* ghost;
    unordered_map<int, Cell*> cells;
    create_ghost(q_neighbor, ghost, cells);

    // Need gamma for this BC
    double gamma = 1.4;

    // Apply boundary condition
    ghost->update_wall(cells, gamma);
    vector<double> ghost_q = ghost->q;
    vector<double> ghost_q_correct = {1, 0, -1, 1e5};

    // Test
    double error = 1e-8;
    for (size_t i = 0; i < ghost_q.size(); i++) {
        ASSERT_NEAR(ghost_q_correct[i], ghost_q[i], error);
    }

}
