#pragma once
#include <cmath>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <cell.h>
#include <face.h>
#include <flux.h>
#include <ghost.h>
#include <inputs.h>
#include <mesh_reader.h>
using std::vector;
using std::unordered_map;
using std::unordered_set;


// Forward declare
class TimeIntegrator;

// Class for storing all cells in the flowfield and for performing high-level
// functions on flowfield cells.
class Flowfield {
    public:
        Flowfield(Inputs, MeshReader);
        void calculate_flux(Flux&);
        void apply_reconstruction();
        void iterate(TimeIntegrator*, Flux&);
        Inputs inputs;
        unordered_map<int, Point> vertices;
        vector<int> connectivity;
        unordered_map<int, Cell*> cells;
        unordered_map<int, Ghost*> ghosts;
        unordered_map<int, Volume*> volumes;
        unordered_set<Face*> faces;
        int n_cells;
        int n_nodes;
        double time;
};
