#include "unit_test.h"

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    ::testing::GTEST_FLAG(filter) = "GlobalTest.bunny";
    int result = RUN_ALL_TESTS();
    return result;
}