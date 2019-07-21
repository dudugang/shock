#include <cmath>
#include <iostream>
#include <vector>
#include <flowfield.h>
using std::cout;
using std::endl;
using std::vector;

Flowfield::Flowfield(Input input) {
    // Store inputs
    this->input = input;

    // Create initial conditions
    for (int i; i < input.n_cells; i++) {

        // Choose right or left of shock tube
        const double pi = 3.1415926535897;
        if (i < 40) {
            cells[i].u = sin(pi*i/40) + 1;
        } else {
            cells[i].u = 1;
        }

        // Find cell center
        cells[i].x = i*input.dx + input.dx/2;

        // Create cell and add to cell map
        cells[i].neighbors = vector<int>{i-1, i+1};
    }

}
/*


        # Add ghost cells
        self.cells[-1] = Ghost(-inputs.dx/2, inputs.u_l, -1, [0])
        self.cells[inputs.n_cells] = Ghost(inputs.n_cells*inputs.dx + inputs.dx/2,
                inputs.u_r, inputs.n_cells, [inputs.n_cells-1])

    def calculate_lambda(self):
        for cell in self.cells.values():
            if cell.type == 'flow':
                cell.lambda_plus  = (cell.u + abs(cell.u))/2
                cell.lambda_minus = (cell.u - abs(cell.u))/2

    def calculate_eigvecs(self):
        for cell in self.cells.values():
            if cell.type == 'flow':
                cell.eigvec_right = 1
                cell.eigvec_left  = 1

    def calculate_flux(self):
        for cell in self.cells.values():
            if cell.type == 'flow':
                # F_i+1/2 = TL^+T^-1 Q_i + TL^-T^-1 Q_i+1
                cell.flux_right = (cell.eigvec_right * cell.lambda_plus  * cell.eigvec_right * cell.u
                                +  cell.eigvec_right * cell.lambda_minus * cell.eigvec_right * self.cells[cell.neighbors[1]].u)
                # F_i-1/2 = TL^+T^-1 Q_i-1 + TL^-T^-1 Q_i
                cell.flux_left  = (cell.eigvec_right * cell.lambda_plus  * cell.eigvec_right * self.cells[cell.neighbors[0]].u
                                +  cell.eigvec_right * cell.lambda_minus * cell.eigvec_right * cell.u)

    def apply_time_integrator(self):
        for cell in self.cells.values():
            if cell.type == 'flow':
                cell.u = cell.u - (self.inputs.dt/self.inputs.dx)*(cell.flux_right - cell.flux_left)

    def update_ghosts(self):
        for cell in self.cells.values():
            if cell.type == 'ghost':
                cell.update(self.cells[cell.neighbors[0]])

# Class representing one cell in the flowfield. A cell stores its volume-
# averaged flowfield data, its geometry, and the cell ID's of its neighbors.
class Cell:
    # Initialize
    def __init__(self, x, u, cell_id, neighbors):
        self.x = x
        self.u = u
        self.cell_id = cell_id
        self.neighbors = neighbors
        self.type = 'flow'



# Class for ghost cells in the boundaries. These are used for enforcing boundary
# conditions.
class Ghost(Cell):
    # Initialize
    def __init__(self, x, u, cell_id, neighbors):
        super().__init__(x, u, cell_id, neighbors)
        self.type = 'ghost'

    # Reflective BC
    def update(self, neighbor):
        self.u = neighbor.u
*/
