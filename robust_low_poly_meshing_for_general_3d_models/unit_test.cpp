#include "unit_test.h"

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    // ::testing::GTEST_FLAG(filter) = "GlobalTest.cache*";
    //::testing::GTEST_FLAG(filter) = "GlobalTest.rotate_test*";
    ::testing::GTEST_FLAG(filter) = "GlobalTest.feature_insertion_cube";
    int result = RUN_ALL_TESTS();
    return result;
}