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
};

#endif