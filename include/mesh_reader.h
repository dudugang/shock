#pragma once
#include <iostream>
#include <string>
#include <H5Cpp.h>
using namespace H5;
using std::cout;
using std::endl;
using std::string;

class MeshReader {
    // TODO: Add destructor for everything in here to prevent memory leaks
    public:
        MeshReader(string);
        void read_mesh();
        string mesh_file;
        double *x_coords;
        double *y_coords;
        int *connectivity;
        int n_cells;
    private:
        template<class T>
        T* read_dataset(H5File, string);
};
