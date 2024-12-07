#ifndef __CHARLES_ROBUST_LOW_MESH_FEATURE_EXTRACTION_H__
#define __CHARLES_ROBUST_LOW_MESH_FEATURE_EXTRACTION_H__

#include "charles_mc33_type.h"

void feature_extraction(
    const std::vector<CustomizedPoint3D>& vertices,
    const std::vector<std::vector<int>>& polygons
);

#endif