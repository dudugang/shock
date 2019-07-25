#include <vector>
#include <unordered_map>
#include <inputs.h>
#include <ghost.h>
#include <flow.h>
#include <face.h>
#include <flux.h>
using std::vector;
using std::unordered_map;

// Class for storing all cells in the flowfield and for performing high-level
// functions on flowfield cells.
class Flowfield {
    public:
        Flowfield(Inputs);
        void calculate_flux();
        void apply_time_integrator();
        Inputs inputs;
        unordered_map<int, Cell*> cells;
        unordered_map<int, Face*> faces;
        double time;
};
