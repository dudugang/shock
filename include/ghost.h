#include "../include/cell.h"
#include <vector>
using std::vector;

class Ghost: public Cell {
    public:
        Ghost(double, double, int, vector<int>);
};
