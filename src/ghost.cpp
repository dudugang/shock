#include <iostream>
#include <ghost.h>
using std::vector;

// Constructor
Ghost::Ghost(double x, vector<double> q, int cell_id, vector<int> neighbor,
    Face *left_face, Face *right_face) : Cell(x, q, cell_id, neighbor,
    left_face, right_face) {
    this->type = "ghost";
}

