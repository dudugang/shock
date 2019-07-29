#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <cell.h>
#include <face.h>
#include <flux.h>
#include <ghost.h>
#include <inputs.h>
using std::vector;
using std::unordered_map;
using std::unordered_set;

// Class for storing all cells in the flowfield and for performing high-level
// functions on flowfield cells.
class Flowfield {
    public:
        Flowfield(Inputs);
        void calculate_flux(Flux&);
        void apply_time_integrator();
        Inputs inputs;
        unordered_set<Volume*> volumes;
        unordered_set<Cell*> cells;
        unordered_set<Ghost*> ghosts;
        unordered_set<Face*> faces;
        unordered_map<int, Volume*> id_to_volume;
        unordered_map<int, Face*> id_to_face;
        double time;
};
