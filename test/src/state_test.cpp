#include <state.h>
#include <gtest/gtest.h>
using std::vector;


TEST(physical_to_conserved_test, ShouldWorkForZeroVelocity) {

    // Create test data
    vector<double> physical = {1, 0, 0, 1e5};
    double gamma = 1.4;

    vector<double> conserved = State::physical_to_conserved(physical, gamma);
    vector<double> conserved_correct = {1, 0, 0, physical[3]/(gamma-1)};

    // Test
    double error = 1e-8;
    for (size_t i = 0; i < physical.size(); i++) {
        ASSERT_NEAR(conserved_correct[i], conserved[i], error) << "Index: " << i;
    }

}


TEST(physical_to_conserved_test, ShouldWorkForNonzeroVelocity) {

    // Create test data
    vector<double> physical = {2, 1, 1, 1e5};
    double gamma = 1.4;

    vector<double> conserved = State::physical_to_conserved(physical, gamma);
    vector<double> conserved_correct = {2, 2, 2, physical[3]/(gamma-1)
        + .5*physical[0]*(physical[1]*physical[1] + physical[2]*physical[2])};

    // Test
    double error = 1e-8;
    for (size_t i = 0; i < physical.size(); i++) {
        ASSERT_NEAR(conserved_correct[i], conserved[i], error) << "Index: " << i;
    }

}


TEST(conserved_to_physical_test, ShouldWorkForZeroVelocity) {

    // Create test data
    vector<double> conserved = {1, 0, 0, 1e5};
    double gamma = 1.4;

    vector<double> physical = State::conserved_to_physical(conserved, gamma);
    vector<double> physical_correct = {1, 0, 0, conserved[3]*(gamma-1)};

    // Test
    double error = 1e-8;
    for (size_t i = 0; i < conserved.size(); i++) {
        ASSERT_NEAR(physical_correct[i], physical[i], error) << "Index: " << i;
    }

}

TEST(conserved_to_physical_test, ShouldWorkForNonzeroVelocity) {

    // Create test data
    vector<double> conserved = {2, 100, 100, 1e5};
    double gamma = 1.4;

    vector<double> physical = State::conserved_to_physical(conserved, gamma);
    vector<double> physical_correct = {2, 50, 50, (gamma-1)*(conserved[3]
        - .5*(conserved[1]*conserved[1] + conserved[2]*conserved[2])/conserved[0])};

    // Test
    double error = 1e-8;
    for (size_t i = 0; i < conserved.size(); i++) {
        ASSERT_NEAR(physical_correct[i], physical[i], error) << "Index: " << i;
    }

}
