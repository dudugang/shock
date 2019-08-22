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
