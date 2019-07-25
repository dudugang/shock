#include <face.h>
using std::cout;
using std::endl;

Face::Face(int face_id, Cell *left_cell, Cell *right_cell) {
    cout << "Creating face with ID " << face_id << endl;
    this->face_id = face_id;
    this->left_cell = left_cell;
    this->right_cell = right_cell;
    flux.resize(3);
}
