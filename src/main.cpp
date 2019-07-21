#include <iostream>
#include <flowfield.h>
using std::cout;
using std::endl;

int main(int argc, char* argv[]) {
    // Get inputs
    Input input;

    // Initialize flowfield
    Flowfield flow(input);

    return 0;
}
