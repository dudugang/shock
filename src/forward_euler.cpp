#include <forward_euler.h>
#include <cell.h>
#include <face.h>


// Constructor
ForwardEuler::ForwardEuler() {

    flux_integral.resize(4);

}


// Use the Forward Euler method to integrate in time
void ForwardEuler::integrate(Flowfield &flow, unordered_map<int, Cell*> &cells, Flux& flux, double dt) {

    // Reconstruct flowfield
    flow.apply_reconstruction();

    // Calculate fluxes
    flow.calculate_flux(flux);

    // Loop over every pair of cell IDs and cells
    for (auto &pair : cells) {

        // Define convenient pointer
        Cell* cell = pair.second;

        // Coefficient in front of flux integral
        double coefficient = -(dt/cell->volume);

        // Apply to every equation
        for (int i = 0; i < 4; i++) {

            // Integrate flux over every cell face
            for (auto &face : cell->faces) {

                // If normal vector of face points away from cell, then flux is
                // positive
                if (cell->id == face->neighbors[0]) {
                    flux_integral[i] += face->flux[i] * face->area;
                } else {
                    flux_integral[i] -= face->flux[i] * face->area;
                }

            }

            // Integrate in time
            cell->q[i] += coefficient * flux_integral[i];

            }

    }

}
