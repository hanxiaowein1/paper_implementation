#ifndef __ROBUST_LOW_POLY_MESHING_COMMON_H__
#define __ROBUST_LOW_POLY_MESHING_COMMON_H__

#include <vector>

#include <Eigen/Dense>
std::vector<std::vector<int>> convert_triangles_to_polygons(const std::vector<Eigen::Vector3i>& triangles);


#endif