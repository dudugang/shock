#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <point.h>
using std::cout;
using std::endl;
using std::string;
using std::vector;

// Forward declare
class Face;

// Class for one cell in the flowfield. Stores geometric data and flowfield
// variables, as well as its cell ID and the cell IDs of its neighbors.
class Volume {
    public:
        Volume(vector<Point>, vector<double>, vector<Face*>, int);
        Point center;
        vector<Point> vertices;
        vector<double> q;
        vector<Face*> faces;
        int id;
        double volume;
        string type;
};
