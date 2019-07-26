#include <fstream>
#include <iostream>
#include <flowfield.h>


// Class for outputting information, both to the terminal during runs and also
// to a solution file for postprocessing.
class Output {
    public:
        Output(Inputs);
        void print(Flowfield, int);
        void final_print(Flowfield);
        void write(Flowfield);
};
