#include <vector>
#include <cell.h>
#include <face.h>
using std::vector;

// Inherits class Cell, similar to a regular Cell but has only one neighbor and
// is not included during either space or time integration.
class Ghost: public Cell {
    public:
        Ghost(double, vector<double>, int, vector<Cell*>, Face*, Face*);
        void update();
};
