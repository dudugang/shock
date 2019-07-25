#pragma once
#include <iostream>
#include <vector>
#include <cell.h>
using std::vector;

class Face {
    public:
        Face(int, Cell*, Cell*);
        int face_id;
        Cell *left_cell;
        Cell *right_cell;
        vector<double> flux;
};
