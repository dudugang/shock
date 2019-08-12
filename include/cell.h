#pragma once
#include <vector>
#include <geometry.h>
#include <point.h>
#include <volume.h>
using std::vector;

// Forward declare
class Face;

// Inherits class Volume, this is a flowfield cell that has two neighbors and is
// used during either space and time integration.
class Cell : public Volume {
    public:
        Cell(vector<Point>, vector<double>, vector<Face*>, int);
        Point center;
};
