#pragma once
#include <utility>
#include <iostream>
#include <map>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <H5Cpp.h>
#include <inputs.h>
#include <point.h>
using namespace H5;
using std::cout;
using std::endl;
using std::map;
using std::pair;
using std::string;
using std::unordered_map;
using std::unordered_set;

// Forward declare
class Cell;
class Ghost;
class Face;

class MeshReader {
    // TODO: Add destructor for everything in here to prevent memory leaks
    public:
        MeshReader(string);
        void create_mesh(Inputs&);
        string mesh_file;
        int n_nodes;
        int n_cells;
        unordered_map<int, Point> vertices;
        unordered_map<int, Cell*> cells;
        unordered_map<int, Ghost*> ghosts;
        unordered_set<Face*> faces;

    private:
        void read_hdf5();
        template<class T>
        T* read_dataset(H5File, string);
        void read_dataset_contents(int*, DataSet, DataSpace, DataSpace);
        void read_dataset_contents(double*, DataSet, DataSpace, DataSpace);
        double *x_coords;
        double *y_coords;
        int *connectivity;
        unordered_map<string, int*> bc_connectivity;
        unordered_map<string, int> bc_face_count;
};
