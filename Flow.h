#include <array>
#include <vector>
using std::array;
using std::vector;

class Flow {

    private:
        vector<double> grid;
        vector<double> u;
        array<double,2> boundary_u;
        int n_points;
        double dt;
        double time;
        int n_iter;
        double length;
        double s;
        double alpha;
        bool implicit;

    public:
        Flow();
        void initialize();
        void solve();
        void iterateExplicit(vector<double>&, double&, double, double, double);
        void iterateImplicit(vector<double>&, double&, double, double, double);
        vector<double> solveLinearSystem(vector<vector<double>>&, vector<double>&);
        vector<double> multiplyMatrices(vector<vector<double>>&, vector<double>&);
        vector<double> subtractMatrices(vector<double>&, vector<double>&);
        void output();

};
