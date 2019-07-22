#include <ghost.h>
#include <vector>
using std::vector;

// Constructor
Ghost::Ghost(double x, vector<double> q, int cell_id, vector<int> neighbor)
    : Cell(x, q, cell_id, neighbor) {
    this->type = "flow";
}

