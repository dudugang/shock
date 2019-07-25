#include <vector>
#include <cell.h>
#include <face.h>
using std::vector;

// Inherits class Cell, this is a flowfield cell that has two neighbors and is
// used during either space and time integration.
class Flow: public Cell {
    public:
        Flow(double, vector<double>, int, vector<Cell*>, Face*, Face*);
        void update();
};
