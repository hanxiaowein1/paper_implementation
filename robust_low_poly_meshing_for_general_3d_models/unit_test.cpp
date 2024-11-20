#include "unit_test.h"

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    ::testing::GTEST_FLAG(filter) = "GlobalTest.feature_point_triangle_same_length";
    int result = RUN_ALL_TESTS();
    return result;
}