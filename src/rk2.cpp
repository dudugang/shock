#include <rk2.h>
#include <cell.h>
#include <face.h>


// Constructor
RK2::RK2() {

    flux_integral.resize(4);

}


// Use the 2-stage Runge-Kutta method to integrate in time
void RK2::integrate(Flowfield &flow, unordered_map<int, Cell*> &cells, Flux& flux, double dt) {

    // Vector to store k1
    vector<vector<double> > k1;
    k1.resize(flow.n_cells);
    for (int i = 0; i < flow.n_cells; i++) {
        k1[i] = {0, 0, 0, 0};
    }
    // Vector to store k2
    vector<vector<double> > k2;
    k2.resize(flow.n_cells);
    for (int i = 0; i < flow.n_cells; i++) {
        k2[i] = {0, 0, 0, 0};
    }

    /* -------------------- STAGE 1 -------------------- */

    // Reconstruct flowfield
    flow.apply_reconstruction();

    // Calculate fluxes
    flow.calculate_flux(flux);

    // Loop over every pair of cell IDs and cells
    for (auto &pair : cells) {

        // Define convenient pointer
        Cell* cell = pair.second;

        // Coefficient in front of flux integral
        double coefficient = -1/cell->volume;

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

            // Calculate k1
            k1[cell->id-1][i] = coefficient * flux_integral[i];

            // Update cell center value
            cell->q[i] += k1[cell->id-1][i] * dt;

            }

    }

    /* -------------------- STAGE 2 -------------------- */

    // Reconstruct flowfield
    flow.apply_reconstruction();

    // Calculate fluxes
    flow.calculate_flux(flux);

    // Loop over every pair of cell IDs and cells
    for (auto &pair : cells) {

        // Define convenient pointer
        Cell* cell = pair.second;

        // Coefficient in front of flux integral
        double coefficient = -1/cell->volume;

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

            // Calculate k1
            k2[cell->id-1][i] = coefficient * flux_integral[i];

            // Update cell center value
            cell->q[i] -= k1[cell->id-1][i] * dt;
            cell->q[i] += (k1[cell->id-1][i]/2 + k2[cell->id-1][i]/2) * dt;

            }

    }

}
