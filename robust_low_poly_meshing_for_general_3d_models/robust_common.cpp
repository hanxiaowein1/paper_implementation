#include "robust_common.h"


std::vector<std::vector<int>> convert_triangles_to_polygons(const std::vector<Eigen::Vector3i>& triangles)
{
    std::vector<std::vector<int>> polygons;
    for(const auto& triangle: triangles)
    {
        std::vector<int> polygon(3, 0);
        polygon[0] = triangle[0];
        polygon[1] = triangle[1];
        polygon[2] = triangle[2];
        polygons.emplace_back(std::move(polygon));
    }
    return polygons;
}