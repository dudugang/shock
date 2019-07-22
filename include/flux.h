#include <vector>
using std::vector;

// Class containing static methods for computing fluxes between cells.
class Flux {
    public:
        static vector<double> steger_warming(vector<double>, vector<double>, double);
};
