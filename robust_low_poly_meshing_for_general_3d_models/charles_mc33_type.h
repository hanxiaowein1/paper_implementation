#ifndef __CHARLES_ROBUST_LOW_MESH_TYPES_H__
#define __CHARLES_ROBUST_LOW_MESH_TYPES_H__
#include "charles_mesh.h"

enum class Vertex {v0, v1, v2, v3, v4, v5, v6, v7, vc};

enum class Edge {e0, e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11};

enum class Face {f0, f1, f2, f3, f4, f5};

enum class Axis {x, y, z};

class CustomizedPoint3D : public charles_mesh::Point3D
{
public:
    bool is_feature_point = false;
    double dot(const CustomizedPoint3D& point) const
    {
        return this->x * point.x + this->y * point.y + this->z * point.z;
    }
    CustomizedPoint3D cross(const CustomizedPoint3D& point) const
    {
        CustomizedPoint3D result;
        result.x = this->y * point.z - this->z * point.y;
        result.y = this->z * point.x - this->x * point.z;
        result.z = this->x * point.y - this->y * point.x;
        return result;
    }
    CustomizedPoint3D operator-(const CustomizedPoint3D& point) const
    {
        CustomizedPoint3D result;
        result.x = this->x - point.x;
        result.y = this->y - point.y;
        result.z = this->z - point.z;
        return result;
    }
    CustomizedPoint3D operator+(const CustomizedPoint3D& point) const
    {
        CustomizedPoint3D result;
        result.x = this->x + point.x;
        result.y = this->y + point.y;
        result.z = this->z + point.z;
        return result;
    }
    CustomizedPoint3D operator*(double value) const
    {
        CustomizedPoint3D result;
        result.x = this->x * value;
        result.y = this->y * value;
        result.z = this->z * value;
        return result;
    }
    CustomizedPoint3D operator/(double value) const
    {
        if(value == 0.0f)
        {
            throw std::exception("value cannot be zero!");
        }
        CustomizedPoint3D result;
        result.x = this->x / value;
        result.y = this->y / value;
        result.z = this-> z / value;
        return result;
    }
    bool between(const CustomizedPoint3D& point1, const CustomizedPoint3D& point2)
    {
        if(!(charles_math::in_interval(this->x + charles_mesh::deviation, point1.x, point2.x) && charles_math::in_interval(this->x - charles_mesh::deviation, point1.x, point2.x)))
        {
            return false;
        }
        if(!(charles_math::in_interval(this->y + charles_mesh::deviation, point1.y, point2.y) && charles_math::in_interval(this->y - charles_mesh::deviation, point1.y, point2.y)))
        {
            return false;
        }
        if(!(charles_math::in_interval(this->z + charles_mesh::deviation, point1.z, point2.z) && charles_math::in_interval(this->z - charles_mesh::deviation, point1.z, point2.z)))
        {
            return false;
        }
        return true;
    }
    bool operator==(const CustomizedPoint3D& point) const
    {
        if(!(this->x == point.x))
        {
            return false;
        }
        if(!(this->y == point.y))
        {
            return false;
        }
        if(!(this->z == point.z))
        {
            return false;
        }
        return true;
    }
};

namespace std {
    template<>
    struct hash<CustomizedPoint3D> {
        size_t operator()(const CustomizedPoint3D& point) const {
            size_t x_hash = std::hash<double>()(point.x);
            size_t y_hash = std::hash<double>()(point.y) << 1;
            size_t z_hash = std::hash<double>()(point.z) << 2;
            return x_hash ^ y_hash ^ z_hash;
        }
    };
    template <> 
    struct hash<std::pair<CustomizedPoint3D, CustomizedPoint3D>> {
        inline size_t operator()(const std::pair<CustomizedPoint3D, CustomizedPoint3D>& v) const {
            std::hash<CustomizedPoint3D> point_hasher;
            return point_hasher(v.first) ^ point_hasher(v.second);
        }
    };
}

#endif