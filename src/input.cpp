#include <input.h>

Input::Input() {
    dt = .0001;
    n_cells = 100;
    n_iterations = 5000;
    length = 1;
    u_l = 1;
    u_r = .5;
    dx = length/n_cells;
}
