#include <memory>

#include "feature_extraction.h"
#include "charles_mesh.h"
#include "charles_mc33_type.h"
#include "unit_test.h"
#include "mesh_factory.h"

#include <Eigen/Dense>

/**
 * @brief check if edge will connect feature point after flip
 * this following is for reference, assume that e2 passed in
 *     v1
 *    /  ^
 * e1/    \e3
 *  v      \
 * v2  -e2-> v3
 *  \ <-e4- ^
 * e5\     / e6
 *    v   /
 *     v4
 * @param half_edge 
 * @return true 
 * @return false 
 */
bool feature_point_connected_after_flip(std::shared_ptr<charles_mesh::HalfEdge<CustomizedPoint3D>> half_edge)
{
    auto e2 = half_edge;
    auto e4 = half_edge->opposite;
    if(!e2->next->vertex->position.is_feature_point)
    {
        return false;
    }
    if(!e4->next->vertex->position.is_feature_point)
    {
        return false;
    }
    return true;
}

/**
 * @brief feature extraction algorithm in paper
 * this following is for reference
 *     v1
 *    /  ^
 * e1/    \e3
 *  v      \
 * v2  -e2-> v3
 *  \ <-e4- ^
 * e5\     / e6
 *    v   /
 *     v4
 * 
 * @param vertices 
 * @param polygons 
 */
void feature_extraction(
    const std::vector<CustomizedPoint3D>& vertices,
    const std::vector<std::vector<int>>& polygons
)
{
    // convert to half edge data structure
    std::shared_ptr<charles_mesh::Mesh<CustomizedPoint3D>> mesh(new charles_mesh::Mesh<CustomizedPoint3D>(vertices, polygons));
    // edge flip between feature point, and avoid self intersection
    for(auto half_edge: mesh->half_edges)
    {
        if(feature_point_connected_after_flip(half_edge))
        {
            // do edge flip that avoid inner interpolation
            mesh->edge_flip_with_intersection_detect(half_edge);
        }
    }
    // save it to obj file to check the result(just for test)
    std::string bunny_obj_src_path = "./bunny.obj";
    charles_mesh::ObjMeshIO<CustomizedPoint3D> obj_mesh_io;
    // auto mesh = obj_mesh_io.load_mesh(bunny_obj_src_path);
    obj_mesh_io.save_mesh("./", "bunny.obj", mesh);
}


TEST(GlobalTest, feature_extraction)
{
    // feature_extraction();
}