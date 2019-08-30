#include <time_integrator.h>
#include <forward_euler.h>
#include <rk2.h>


// Return an instance of a derived class that is the correct time integrator
// chosen in the inputs
TimeIntegrator* TimeIntegrator::choose_time_integrator(Inputs inputs) {

    if (inputs.time_integrator == "ForwardEuler") {

        ForwardEuler* time_integrator = new ForwardEuler();
        return time_integrator;

    } else if (inputs.time_integrator == "RK2") {

        RK2* time_integrator = new RK2();
        return time_integrator;

    } else {

        cout << "Invalid time integrator chosen!" << endl;

    }

}
