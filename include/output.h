#pragma once
#include <fstream>
#include <iostream>
#include <flowfield.h>
using std::cout;
using std::endl;
using std::ofstream;


// Class for outputting information, both to the terminal during runs and also
// to a solution file for postprocessing.
class Output {
    public:
        Output(Inputs, int);
        void print(Flowfield, int);
        void final_print(Flowfield);
        void write(Flowfield, int);
};
