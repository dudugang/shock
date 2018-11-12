#include <array>
#include <vector>
using std::array;
using std::vector;

class Flow {
    private:
        vector<double> grid;
        vector<vector<double> > q;
        vector<vector<double> > q_vertex;
        vector<vector<double> > f;
        vector<double> p;
        // 3x2 array, with rho, u, and e on the rows and left/right boundaries
        // on the columns
        array<array<double,2>,3> boundary_conditions;
        double gamma;
        int n_cells;
        double cfl;
        double dt;
        double time;
        int n_iter;
        double length;
        double s;

    public:
        Flow();
        void initialize();
        vector<vector<double> > calculate_f_vector(vector<vector<double> >&, int, double);
        vector<double> calculate_pressure(vector<vector<double> >&, double);
        vector<double> calculate_u(vector<vector<double> >&, double);
        double calculate_dt(vector<vector<double> >&, double, double, double);
        void solve();
        void iterate(vector<vector<double> >&, vector<vector<double> >&, double&, double&, double&, double&, int&);
        vector<double> case_1_riemann(double, double, double, double, double, double, double, double);
        vector<double> case_2_riemann(double, double, double, double, double, double, double, double);
        double calculate_pstar(double, double, double, double, double, double, double);
        double f_pstar(double, double, double);
        double f_prime_pstar(double, double, double);
        void write();
        void output();
};
