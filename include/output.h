#include <iostream>
#include <flowfield.h>


// Class for outputting information, both to the terminal during runs and also
// to a solution file for postprocessing.
class Output {
    public:
        void print(Flowfield);
        void write();
};
