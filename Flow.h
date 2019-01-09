#include <Eigen/Dense>
#include <vector>
using namespace Eigen;
using std::vector;

class Flow {
    private:
        vector<double> grid;
        vector<VectorXd> q;
        vector<VectorXd> f_right;
        vector<VectorXd> f_left;
        double gamma;
        double r_gas;
        int n_equations;
        int n_cells;
        int n_ghosts;
        double cfl;
        double max_dt;
        double dt;
        double time;
        int n_iter;
        double length;
        double dx;

    public:
        Flow();
        void initialize();
        //double calculate_dt(vector<vector<double> >&, double, double, double);
        void solve();
        void iterate(vector<VectorXd>&, double&, double, double, double, int);
        void write();
        void output();
};
