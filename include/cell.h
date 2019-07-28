#include <vector>
#include <volume.h>
#include <face.h>
using std::vector;

// Inherits class Volume, this is a flowfield cell that has two neighbors and is
// used during either space and time integration.
class Cell : public Volume {
    public:
        Cell(double, vector<double>, int, vector<Volume*>, Face*, Face*);
        void update();
};
