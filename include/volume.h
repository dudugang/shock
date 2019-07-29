#pragma once
#include <iostream>
#include <string>
#include <vector>
using std::string;
using std::vector;

// Forward declare Face, since Volume and Face depend on each other
class Face;

// Class for one cell in the flowfield. Stores geometric data and flowfield
// variables, as well as its cell ID and the cell IDs of its neighbors.
class Volume {
    public:
        Volume(vector<double>, vector<double>, int, vector<Volume*>, Face*, Face*);
        vector<double> center;
        vector<double> q;
        int volume_id;
        vector<Volume*> neighbors;
        Face *left_face;
        Face *right_face;
        string type;
        virtual void update() = 0;
};
