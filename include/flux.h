#include <vector>
using std::vector;

// Class containing static methods for computing fluxes between cells.
class Flux {
    public:
        Flux();
        vector<double> steger_warming(vector<double>, vector<double>, double);
    private:
        vector<double> eigvalues_l;
        vector<double> eigvalues_r;
        vector<double> eig_split_l;
        vector<double> eig_split_r;
        vector<double> flux_l;
        vector<double> flux_r;
        vector<double> total_flux;
};
