#pragma once
#include <vector>
using std::vector;


// This class defines static functions that give information about the state of
// a gas.
class State {

    public:
        static vector<double> physical_to_conserved(vector<double>, double);

};
