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
    bool between(const CustomizedPoint3D& point1, const CustomizedPoint3D& point2)
    {
        if(!(charles_math::in_interval(this->x, point1.x, point2.x)))
        {
            return false;
        }
        if(!(charles_math::in_interval(this->y, point1.y, point2.y)))
        {
            return false;
        }
        if(!(charles_math::in_interval(this->z, point1.z, point2.z)))
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

#endif