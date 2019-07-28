#pragma once
#include <iostream>
#include <vector>
#include <volume.h>
using std::vector;

class Face {
    public:
        Face(int, Volume*, Volume*);
        int face_id;
        Volume *left_volume;
        Volume *right_volume;
        vector<double> flux;
};
