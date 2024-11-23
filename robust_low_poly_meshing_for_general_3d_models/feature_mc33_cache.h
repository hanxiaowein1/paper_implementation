#ifndef __CHARLES_FEATURE_MC33_CACHE_H__
#define __CHARLES_FEATURE_MC33_CACHE_H__

#include <unordered_map>
#include <boost/functional/hash.hpp>
#include "charles_mc33_type.h"

// key: z, y, x, edge(from 0 ~ 7)
template <typename T>
class FeatureMC33Cache
{
private:
    std::unordered_map<std::tuple<int, int, int, Edge>, T, boost::hash<std::tuple<int, int, int, Edge>>> cache;
public:
    T get(const std::tuple<int, int, int, Edge>& key);
    void set(const std::tuple<int, int, int, Edge>& key, const T& value);
    int size()
    {
        return this->cache.size();
    }
};

template <typename T>
T FeatureMC33Cache<T>::get(const std::tuple<int, int, int, Edge>& key)
{
    // interesting, edge need be recalculated, because (0, 0, 0, Edge::e4) equals (0, 0, 1, Edge::e0)
    const auto& [z, y, x, edge] = key;
    switch(edge)
    {
    case Edge::e0:
        if(x - 1 >= 0 && z - 1 >= 0)
        {
            return this->cache.at({z - 1, y, x - 1, Edge::e6});
        }
        if(z - 1 >= 0)
        {
            return this->cache.at({z - 1, y, x, Edge::e2});
        }
        if(x - 1 >= 0)
        {
            return this->cache.at({z, y, x - 1, Edge::e4});
        }
    case Edge::e1:
        if(x - 1 >= 0)
        {
            return this->cache.at({z, y, x - 1, Edge::e5});
        }
    case Edge::e2:
        if(x - 1 >= 0)
        {
            return this->cache.at({z, y, x - 1, Edge::e6});
        }
    case Edge::e3:
        if(x - 1 >= 0 && y - 1 >= 0)
        {
            return this->cache.at({z, y - 1, x - 1, Edge::e5});
        }
        if(y - 1 >= 0)
        {
            return this->cache.at({z, y - 1, x, Edge::e1});
        }
        if(x - 1 >= 0)
        {
            return this->cache.at({z, y, x - 1, Edge::e7});
        }
    case Edge::e4:
        if(z - 1 >= 0)
        {
            return this->cache.at({z - 1, y, x, Edge::e6});
        }
    case Edge::e5:
        // not possible
        break;
    case Edge::e6:
        // not possible
        break;
    case Edge::e7:
        if(y - 1 >= 0)
        {
            return this->cache.at({z, y - 1, x, Edge::e5});
        }
    }
    return this->cache.at(key);
}

template <typename T>
void FeatureMC33Cache<T>::set(const std::tuple<int, int, int, Edge>& key, const T& value)
{
    this->cache.emplace(key, value);
}

#endif