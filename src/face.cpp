#include <face.h>

Face::Face(int face_id, Cell *left_cell, Cell *right_cell) {
    this->face_id = face_id;
    this->left_cell = left_cell;
    this->right_cell = right_cell;
}
