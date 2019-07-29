#include <vector>
#include <volume.h>
#include <face.h>
using std::vector;

// Inherits class Volume, similar to a regular Cell but has only one neighbor
// and is not included during either space or time integration.
class Ghost : public Volume {
    public:
        Ghost(vector<double>, vector<double>, int, vector<Volume*>, Face*, Face*);
        void update();
};
