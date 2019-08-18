#pragma once
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <H5Cpp.h>
#include <flowfield.h>
using namespace H5;
using std::cout;
using std::endl;
using std::ofstream;
using std::unordered_map;

class Cell;

// Class for outputting information, both to the terminal during runs and also
// to a solution file for postprocessing.
class Output {
    public:
        Output(Inputs, Flowfield&);
        void add_time(double);
        void print(Flowfield, int);
        void final_print(Flowfield);
        void write(Flowfield, int);
        void write_results(string, unordered_map<int, Cell*>&, int);
        void write_dataset(H5File, string, int[], int);
        void write_dataset(H5File, string, double[], int);

    private:
        vector<double> times;
};
