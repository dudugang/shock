#pragma once
#include <iostream>
#include <string>
#include <vector>
using std::string;
using std::vector;

// Forward declare Face, since Cell and Face depend on each other
class Face;

// Class for one cell in the flowfield. Stores geometric data and flowfield
// variables, as well as its cell ID and the cell IDs of its neighbors.
class Cell {
    public:
        Cell(double, vector<double>, int, vector<Cell*>, Face*, Face*);
        double x;
        vector<double> q;
        int cell_id;
        vector<Cell*> neighbors;
        Face *left_face;
        Face *right_face;
        string type;
        virtual void update() = 0;
};
