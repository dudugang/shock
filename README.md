Flow
==========

Flow is a 1D finite-volume CFD code that solves the Euler Equations for
compressible, inviscid flow with energy transfer. A first-order reconstruction
is used within each cell and the fluxes are calculated with the Steger and
Warming Flux Vector Splitting Method. Time integration is done with a forward
Euler time integrator. Options for other numerical methods will be added in the
future. Input values like constants and initial conditions are hard-coded in
Flow.cpp, but an input file will eventually be added. The executable is
'bin/flow', which outputs a file 'solution.dat' containing grid and flow data.
A simple plotting utility 'animation.py' is included for visualizing the
results.

Getting Started
----------

This section will explain how to get a working copy of the code on your machine.

### Prerequisites

The code is written in C++ and compiles with g++. Also, GNU make is used to compile
the code. Many Linux systems come with these installed, but if not, a working
version can be installed with:

    sudo apt-get install g++
    sudo apt-get install make

Other compilers may work (as long as they support C++11) but are not supported
at this time.

The code uses two other libraries, Eigen (for linear algebra) and Google Test
(for unit testing). However, the project includes these libraries for now. Eigen
is header-only and requires no compiling, but Google Test is compiled and linked
with the rest of the code. Currently, a precompiled version of Google Test is
included and should work on other Linux systems, but later on this will be added
to the compile process.

### Installing

To install the code, download the repository from GitHub:

    git clone https://github.com/alilasemi/flow.git

Go to the main directory and compile using the provided build script:

    cd flow
    ./build.sh

To run the program, run the executable:

    bin/flow

This should generate a file solution.dat which contains the results of the run.

Running the Tests
----------

Unit testing is done during the build process at the end of the build script. To
run the tests manually, the test executable can also be run:

    bin/flowtest

All tests currently run successfully, but several portions of the code still
need unit tests to be written.

Author
----------

**Ali Lasemi** - email: alasem2@illinois.edu - GitHub: alilasemi
