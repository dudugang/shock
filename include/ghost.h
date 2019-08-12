#pragma once
#include <unordered_map>
#include <vector>
#include <inputs.h>
#include <point.h>
#include <volume.h>
using std::unordered_map;
using std::vector;

// Forward declare
class Face;
class Cell;

// Inherits class Volume, similar to a regular Cell but has only one neighbor
// and is not included during either space or time integration.
class Ghost : public Volume {
    public:
        Ghost(vector<Point>, vector<double>, vector<Face*>, int);
        void update(Inputs&, unordered_map<int, Cell*>&);
        void update_inflow(unordered_map<int, Cell*>&);
        void update_wall(unordered_map<int, Cell*>&);
};
