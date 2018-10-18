#include <array>
#include <vector>
using std::array;
using std::vector;

class Flow {

    private:
        vector<double> grid;
        vector<vector<double> > q;
        vector<vector<double> > f;
        vector<double> p;
        // 3x2 array, with rho, u, and e on the rows and left/right boundaries
        // on the columns
        array<array<double,2>,3> boundary_conditions;
        double gamma;
        int n_cells;
        double dt;
        double time;
        int n_iter;
        double length;
        double s;

    public:
        Flow();
        void initialize();
        vector<vector<double> > calculate_f_vector(vector<vector<double> >&, int, double);
        vector<double> calculate_pressure(vector<vector<double> >&, int, double);
        void solve();
        void iterate(vector<vector<double> >&, double&, double&, double&, double&, int&);
        void write();
        void output();

};
