/**
 * @file feature_mc33_cache.cpp
 * @author charles
 * @brief the following can be cached
 * 1. an internal edge shared by four cubes, so there are three duplicate computation in interpolation between vertex on same edge. So interpolation on edge can be cached
 * @version 0.1
 * @date 2024-11-22
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include "feature_mc33_cache.h"
#include "unit_test.h"

TEST(GlobalTest, cache_get_e0)
{
    // has diagonal cube
    {
        FeatureMC33Cache<int> feature_mc33_cache;
        feature_mc33_cache.set({0, 0, 0, Edge::e6}, 1);
        auto value = feature_mc33_cache.get({1, 0, 1, Edge::e0});
        ASSERT_EQ(value, 1);
    }

    // has down cube
    {
        FeatureMC33Cache<int> feature_mc33_cache;
        feature_mc33_cache.set({0, 0, 0, Edge::e2}, 1);
        auto value = feature_mc33_cache.get({1, 0, 0, Edge::e0});
        ASSERT_EQ(value, 1);
    }

    // has y left cube
    {
        FeatureMC33Cache<int> feature_mc33_cache;
        feature_mc33_cache.set({0, 0, 0, Edge::e4}, 1);
        auto value = feature_mc33_cache.get({0, 0, 1, Edge::e0});
        ASSERT_EQ(value, 1);
    }

    // throw exception
    {
        EXPECT_THROW({
            FeatureMC33Cache<int> feature_mc33_cache;
            auto value = feature_mc33_cache.get({0, 0, 1, Edge::e0});
        }, std::out_of_range);
    }
}

TEST(GlobalTest, cache_get_e1)
{
    // has x left
    {
        FeatureMC33Cache<int> feature_mc33_cache;
        feature_mc33_cache.set({0, 0, 0, Edge::e5}, 1);
        auto value = feature_mc33_cache.get({0, 0, 1, Edge::e1});
        ASSERT_EQ(value, 1);
    }

    // already have
    {
        FeatureMC33Cache<int> feature_mc33_cache;
        feature_mc33_cache.set({0, 0, 0, Edge::e1}, 1);
        auto value = feature_mc33_cache.get({0, 0, 0, Edge::e1});
        ASSERT_EQ(value, 1);
    }
}

TEST(GlobalTest, cache_get_e2)
{
    // has x left
    {
        FeatureMC33Cache<int> feature_mc33_cache;
        feature_mc33_cache.set({0, 0, 0, Edge::e6}, 1);
        auto value = feature_mc33_cache.get({0, 0, 1, Edge::e2});
        ASSERT_EQ(value, 1);
    }
}

TEST(GlobalTest, cache_get_e3)
{
    // has diagonal cube
    {
        FeatureMC33Cache<int> feature_mc33_cache;
        feature_mc33_cache.set({0, 0, 0, Edge::e5}, 1);
        auto value = feature_mc33_cache.get({0, 1, 1, Edge::e3});
        ASSERT_EQ(value, 1);
    }

    // has y left
    {
        FeatureMC33Cache<int> feature_mc33_cache;
        feature_mc33_cache.set({0, 0, 0, Edge::e1}, 1);
        auto value = feature_mc33_cache.get({0, 1, 0, Edge::e3});
        ASSERT_EQ(value, 1);
    }

    // has x left
    {
        FeatureMC33Cache<int> feature_mc33_cache;
        feature_mc33_cache.set({0, 0, 0, Edge::e7}, 1);
        auto value = feature_mc33_cache.get({0, 0, 1, Edge::e3});
        ASSERT_EQ(value, 1);
    }
}

TEST(GlobalTest, cache_get_e4)
{
    // has downside cube
    {
        FeatureMC33Cache<int> feature_mc33_cache;
        feature_mc33_cache.set({0, 0, 0, Edge::e6}, 1);
        auto value = feature_mc33_cache.get({1, 0, 0, Edge::e4});
        ASSERT_EQ(value, 1);
    }
}

TEST(GlobalTest, cache_get_e7)
{
    // has y left cube
    {
        FeatureMC33Cache<int> feature_mc33_cache;
        feature_mc33_cache.set({0, 0, 0, Edge::e5}, 1);
        auto value = feature_mc33_cache.get({0, 1, 0, Edge::e7});
        ASSERT_EQ(value, 1);
    }
}