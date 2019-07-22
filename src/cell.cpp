#include <cell.h>
#include <vector>
using std::vector;

// Constructor
Cell::Cell(double x, vector<double> q, int cell_id, vector<int> neighbors) {
    this->x = x;
    this->q = q;
    this->cell_id = cell_id;
    this->neighbors = neighbors;
    this->type = "flow";
}
