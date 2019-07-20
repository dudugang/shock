#include "../include/ghost.h"
#include <vector>
using std::vector;

// Constructor
Ghost::Ghost(double x, double u, int cell_id, vector<int> neighbors) {
    this->type = "flow";
}
