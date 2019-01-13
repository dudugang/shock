#include <gtest/gtest.h>

int main(int argc, char **argv) {
    // Run all tests
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
