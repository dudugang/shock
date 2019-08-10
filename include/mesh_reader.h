#pragma once
#include <utility>
#include <iostream>
#include <map>
#include <string>
#include <unordered_map>
#include <H5Cpp.h>
#include <cell.h>
#include <point.h>
using namespace H5;
using std::cout;
using std::endl;
using std::map;
using std::pair;
using std::string;
using std::unordered_map;

class MeshReader {
    // TODO: Add destructor for everything in here to prevent memory leaks
    public:
        MeshReader(string);
        void create_mesh();
        string mesh_file;
        int n_nodes;
        int n_cells;
        unordered_map<int, Point> vertices;
        unordered_map<int, Cell*> cells;

    private:
        void read_hdf5();
        template<class T>
        T* read_dataset(H5File, string);
        double *x_coords;
        double *y_coords;
        int *connectivity;
};
