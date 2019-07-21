#include <string>
#include <vector>
using std::string;
using std::vector;

class Cell {
    public:
        double x;
        double u;
        int cell_id;
        vector<int> neighbors;
        string type;
};
