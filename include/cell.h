#pragma once
#include <string>
#include <vector>
using std::string;
using std::vector;

// Class for one cell in the flowfield. Stores geometric data and flowfield
// variables, as well as its cell ID and the cell IDs of its neighbors.
class Cell {
    public:
        Cell(double, vector<double>, int, vector<int>);
        double x;
        vector<double> q;
        int cell_id;
        vector<int> neighbors;
        string type;
};
