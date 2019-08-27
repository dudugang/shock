#include <state.h>


// Takes a vector of physical variables and converts them to a vector of
// conserved variables.
vector<double> State::physical_to_conserved(vector<double> physical, double gamma) {

    vector<double> conserved = {
        physical[0],
        physical[0] * physical[1],
        physical[0] * physical[2],
        physical[3]/(gamma-1) + .5*physical[0]*(physical[1]*physical[1]
                                              + physical[2]*physical[2])
    };

    return conserved;

}


// Takes a vector of conserved variables and converts them to a vector of
// physical variables.
vector<double> State::conserved_to_physical(vector<double> conserved, double gamma) {

    vector<double> physical = {
        conserved[0],
        conserved[1] / conserved[0],
        conserved[2] / conserved[0],
       (conserved[3] - .5*(conserved[1]*conserved[1]
            + conserved[2]*conserved[2])/conserved[0]) * (gamma-1)
    };

    return physical;

}
