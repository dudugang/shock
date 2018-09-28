#include <array>
#include <vector>
using std::array;
using std::vector;

class Flow {

    private:
        vector<double> grid;
        vector<vector<double>,3> q;
        // 3x2 array, with rho, u, and e on the rows and left/right boundaries
        // on the columns
        array<array<double,2>,3> boundary_conditions;
        int n_cells;
        double dt;
        double time;
        int n_iter;
        double length;
        double s;

    public:
        Flow();
        void initialize();
        void solve();
        void iterateExplicit(vector<double>&, double&, double, double, double);
        void iterateImplicit(vector<double>&, double&, double, double, double);
        vector<double> solveLinearSystem(vector<vector<double>>&, vector<double>&);
        vector<double> multiplyMatrices(vector<vector<double>>&, vector<double>&);
        vector<double> subtractMatrices(vector<double>&, vector<double>&);
        void write();
        void output();

};
