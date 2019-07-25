#include <cell.h>
#include <vector>
using std::cout;
using std::endl;
using std::vector;

// Constructor
Cell::Cell(double x, vector<double> q, int cell_id, vector<Cell*> neighbors,
    Face *left_face, Face *right_face) {
    cout << "Creating cell with ID " << cell_id << endl;
    this->x = x;
    this->q = q;
    this->cell_id = cell_id;
    this->neighbors = neighbors;
    this->left_face = left_face;
    this->right_face = right_face;
}
