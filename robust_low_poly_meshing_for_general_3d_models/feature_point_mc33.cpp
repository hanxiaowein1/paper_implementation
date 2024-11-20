#include <vector>
#include <iterator>
#include <numeric>
#include <fstream>
#include <unordered_set>
#include <bitset>

#include "Eigen/Dense"
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include "igl/readOFF.h"
#include "igl/signed_distance.h"

#include "feature_point_mc33.h"
#include "feature_point_mc33_tables.h"
#include "unit_test.h"
#include "charles_mc33_type.h"


double interpolate_between_four_edges(
    const std::vector<Eigen::Vector3d>& coors,
    const std::vector<double>& signed_distance,
    const Axis& axis
)
{
    std::vector<Edge> edges;
    std::vector<double> edges_interpolation(4, 0.0f);
    double interpolation_axis = 0.0f;
    switch(axis)
    {
    case Axis::x:
        edges = std::vector<Edge>{Edge::e8, Edge::e9, Edge::e10, Edge::e11};
        break;
    case Axis::y:
        edges = std::vector<Edge>{Edge::e0, Edge::e2, Edge::e4, Edge::e6};
        break;
    case Axis::z:
        edges = std::vector<Edge>{Edge::e1, Edge::e3, Edge::e5, Edge::e7};
        break;
    default:
        throw std::invalid_argument("invalid argument!");
    }
    for(auto iterator = edges.begin(); iterator != edges.end(); iterator++)
    {
        int index = std::distance(edges.begin(), iterator);
        auto edge = *iterator;
        // get vertex from edge
        auto vertices = EDGE_VERTEX.at(edge);
        auto first_point = *(vertices.begin());
        auto second_point = *(vertices.begin()++);
        double t = signed_distance[static_cast<int>(first_point)] / (signed_distance[static_cast<int>(first_point)] - signed_distance[static_cast<int>(second_point)]);
        edges_interpolation[index] = coors[static_cast<int>(first_point)].x() + t * std::abs(coors[static_cast<int>(first_point)].x() - coors[static_cast<int>(second_point)].x());
    }
    interpolation_axis = std::accumulate(edges_interpolation.begin(), edges_interpolation.end(), 0.0f);

    return interpolation_axis;
}

// get interpolation point by all points on cube
Eigen::Vector3d get_vertex(
    const std::vector<Eigen::Vector3d>& coors,
    const std::vector<double>& signed_distance
)
{
    // interpolate between x(edge 8 9 A B)
    double interpolation_x = interpolate_between_four_edges(coors, signed_distance, Axis::x);
    double interpolation_y = interpolate_between_four_edges(coors, signed_distance, Axis::y);
    double interpolation_z = interpolate_between_four_edges(coors, signed_distance, Axis::z);
    Eigen::Vector3d vertex;
    vertex[0] = interpolation_x;
    vertex[1] = interpolation_y;
    vertex[2] = interpolation_z;
    return vertex;
}

// interpolation between two points on same edge
Eigen::Vector3d get_vertex(
    const std::vector<Eigen::Vector3d>& coors,
    const std::vector<double>& signed_distance,
    bool on_x, bool on_y, bool on_z,
    Vertex first_point, Vertex second_point
)
{
    Eigen::Vector3d vertex;
    double t = signed_distance[static_cast<int>(first_point)] / (signed_distance[static_cast<int>(first_point)] - signed_distance[static_cast<int>(second_point)]);
    if(on_x)
    {
        vertex[0] = coors[static_cast<int>(first_point)].x() + t * std::abs(coors[static_cast<int>(first_point)].x() - coors[static_cast<int>(second_point)].x());
    }
    else
    {
        vertex[0] = coors[static_cast<int>(first_point)].x();
    }
    if(on_y)
    {
        vertex[1] = coors[static_cast<int>(first_point)].y() + t * std::abs(coors[static_cast<int>(first_point)].y() - coors[static_cast<int>(second_point)].y());
    }
    else
    {
        vertex[1] = coors[static_cast<int>(first_point)].y();
    }
    if(on_z)
    {
        vertex[2] = coors[static_cast<int>(first_point)].z() + t * std::abs(coors[static_cast<int>(first_point)].z() - coors[static_cast<int>(second_point)].z());
    }
    else
    {
        vertex[2] = coors[static_cast<int>(first_point)].z();
    }
    return vertex;
}

Eigen::Vector3d get_vertex_by_edge(
    Edge edge,
    const std::vector<Eigen::Vector3d>& coors,
    const std::vector<double>& signed_distance
)
{
    Eigen::Vector3d vertex;
    switch (edge)
    {
    case Edge::e0:
        // vertex 0 1
        vertex = get_vertex(coors, signed_distance, false, true, false, Vertex::v0, Vertex::v1);
        break;
    case Edge::e1:
        // vertex 1 2
        vertex = get_vertex(coors, signed_distance, false, false, true, Vertex::v1, Vertex::v2);
        break;
    case Edge::e2:
        // vertex 2 3
        vertex = get_vertex(coors, signed_distance, false, true, false, Vertex::v3, Vertex::v2);
        break;
    case Edge::e3:
        // vertex 0 3
        vertex = get_vertex(coors, signed_distance, false, false, true, Vertex::v0, Vertex::v3);
        break;
    case Edge::e4:
        // vertex 4 5
        vertex = get_vertex(coors, signed_distance, false, true, false, Vertex::v4, Vertex::v5);
        break;
    case Edge::e5:
        // vertex 5 6
        vertex = get_vertex(coors, signed_distance, false, false, true, Vertex::v5, Vertex::v6);
        break;
    case Edge::e6:
        // vertex 7 6
        vertex = get_vertex(coors, signed_distance, false, true, false, Vertex::v7, Vertex::v6);
        break;
    case Edge::e7:
        // vertex 4 7
        vertex = get_vertex(coors, signed_distance, false, false, true, Vertex::v4, Vertex::v7);
        break;
    case Edge::e8:
        // vertex 0 4
        vertex = get_vertex(coors, signed_distance, true, false, false, Vertex::v0, Vertex::v4);
        break;
    case Edge::e9:
        // vertex 1 5
        vertex = get_vertex(coors, signed_distance, true, false, false, Vertex::v1, Vertex::v5);
        break;
    case Edge::e10:
        // vertex 2 6
        vertex = get_vertex(coors, signed_distance, true, false, false, Vertex::v2, Vertex::v6);
        break;
    case Edge::e11:
        // vertex 3 7
        vertex = get_vertex(coors, signed_distance, true, false, false, Vertex::v3, Vertex::v7);
        break;
    //case Edge::ec:
    //    // center interpolated point
    //    vertex = get_vertex(coors, signed_distance);
    default:
        throw std::invalid_argument("invalid edge!");
        // break;
    }
    return vertex;
}

// std::vector<unsigned short> get_edgess(const std::vector<double>& signed_distances)
// {
//     std::bitset<8> distances_signs{0b11111111};
//     for(int i = 0; i < signed_distances.size(); i++)
//     {
//         auto signed_distance = signed_distances[i];
//         if(signed_distance < 0)
//         {
//             distances_signs.set(i, false);
//         }
//     }
//     auto sub_cases = MC33_TABLES.at(distances_signs);
//     // get all connected vertices that are not ambiguous
//     std::unordered_set<std::unordered_set<Vertex>, boost::hash<std::unordered_set<Vertex>>> connected_verticess;
//     for(unsigned short i = 0; i < 8; i++)
//     {
//         for(unsigned short j = i + 1; j < 8; j++)
//         {
//             Vertex vertex1 = static_cast<Vertex>(i), vertex2 = static_cast<Vertex>(j);
//             if (distances_signs[static_cast<unsigned short>(vertex1)] == distances_signs[static_cast<unsigned short>(vertex2)])
//             {
//                 if (!vertex_connected(distances_signs, vertex1, vertex2))
//                 {
//                     if (vertex_interpolation_connected(signed_distances, vertex1, vertex2))
//                     {
//                         connected_verticess.emplace(std::unordered_set<Vertex>{vertex1, vertex2});
//                     }
//                 }
//             }
//         }
//     }

//     // get egess by connected vertices
//     auto triangles = sub_cases.at(connected_verticess);
//     return triangles;
// }

void write_obj(std::string filename, const std::vector<Eigen::Vector3d>& vertices, const std::vector<Eigen::Vector3i>& triangles)
{
	std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file for writing: " << filename << std::endl;
        return;
    }
    for(const auto& vertice: vertices)
    {
        file << "v " << vertice[0] << " " << vertice[1] << " " << vertice[2] << std::endl;
    }
    for(const auto& triangle: triangles)
    {
		file << "f " << triangle[0] + 1 << " " << triangle[1] + 1 << " " << triangle[2] + 1 << std::endl;
    }
	file.close();
}

Eigen::Vector3d get_interpolated_vertex_normal(const Eigen::Vector3d& interpolated_vertex)
{
    // TODO: need implementation
    Eigen::Vector3d normal;
    return normal;
}

Eigen::Vector3d get_feature_vertex(std::vector<Eigen::Vector3d>& interpolated_vertices, std::vector<Eigen::Vector3d>& normals)
{
    // TODO: need implamentation, should use multivariable minimum/maximum function
    Eigen::Vector3d feature_vertex;
    return feature_vertex;
}

Eigen::Vector3d get_feature_vertex(
    const std::tuple<unsigned short, unsigned short>& feature_point,
    const std::vector<Eigen::Vector3d>& coors,
    const std::vector<double>& signed_distances
)
{
    auto [length, edges] = feature_point;
    std::vector<Eigen::Vector3d> interpolated_vertices;
    std::vector<Eigen::Vector3d> normals;
    for(int i = 0; i < length; i++)
    {
        edges = edges >> (i == 0 ? 0 : 4);
        Edge edge = static_cast<Edge>(edges & 0xF);
        auto vertex = get_vertex_by_edge(edge, coors, signed_distances);
        auto normal = get_interpolated_vertex_normal(vertex);
        interpolated_vertices.emplace_back(std::move(vertex));
        normals.emplace_back(normal);
    }
    // get feature point from normals and interpolated vertices
    auto feature_vertex = get_feature_vertex(interpolated_vertices, normals);
    return feature_vertex;
}


std::tuple<std::vector<Eigen::Vector3d>, std::vector<Eigen::Vector3i>> get_feature_mc33_triangles(
    const FeatureMC33Table& feature_mc33_table,
    const std::vector<Eigen::Vector3d>& coors,
    const std::vector<double>& signed_distances
)
{
    std::vector<Eigen::Vector3d> vertices;
    std::vector<Eigen::Vector3i> triangles;

    std::vector<Eigen::Vector3d> feature_vertices;
    for(const auto& feature_point: feature_mc33_table.feature_points)
    {
        auto feature_vertice = get_feature_vertex(feature_point, coors, signed_distances);
        feature_vertices.emplace_back(feature_vertice);
    }
    auto lambda_func = [&](const Edge& edge) -> Eigen::Vector3d {
        auto vertex = get_vertex_by_edge(edge, coors, signed_distances);
        return vertex;
    };
    std::function<Eigen::Vector3d(const Edge&)> func = lambda_func;
    // compute feature triangles
    for(const auto& [feature_index, f_triangles]: feature_mc33_table.feature_triangles)
    {
        for(const auto& f_triangle: f_triangles)
        {
            auto temp_vertices = handle_edges(f_triangle, FeatureMC33Table::feature_triangle_length, func);
            temp_vertices.emplace_back(feature_vertices[feature_index]);
            if(temp_vertices.size() != 3)
            {
                throw std::exception("triangle that contains feature vertex doesn't has three vertice");
            }
            int current_vertex_index = vertices.size() - 1;
            triangles.emplace_back(Eigen::Vector3i{current_vertex_index, current_vertex_index + 1, current_vertex_index + 2});
            vertices.insert(vertices.end(), temp_vertices.begin(), temp_vertices.end());
        }
    }

    // compute with mc33 inserted point triangles
    // get center interpolated point
    auto mc33_point = get_vertex(coors, signed_distances);
    for(const auto& mc33_triangle: feature_mc33_table.mc33_triangles)
    {
        // auto lambda_func = [&](const Edge& edge) -> Eigen::Vector3d
        auto temp_vertices = handle_edges(mc33_triangle, FeatureMC33Table::mc33_triangle_length, func);
        temp_vertices.emplace_back(mc33_point);
        if(temp_vertices.size() != 3)
        {
            throw std::exception("triangle that contains mc33 point doesn't has three vertice");
        }
        int current_vertex_index = vertices.size() - 1;
        triangles.emplace_back(Eigen::Vector3i{current_vertex_index, current_vertex_index + 1, current_vertex_index + 2});
        vertices.insert(vertices.end(), temp_vertices.begin(), temp_vertices.end());
    }

    // compute normal triangles
    for(const auto& common_triangle: feature_mc33_table.common_triangles)
    {
        auto temp_vertices = handle_edges(common_triangle, FeatureMC33Table::common_triangle_length, func);
        if(temp_vertices.size() != 3)
        {
            throw std::exception("common triangle doesn't has three vertice");
        }
        int current_vertex_index = vertices.size() - 1;
        triangles.emplace_back(Eigen::Vector3i{current_vertex_index, current_vertex_index + 1, current_vertex_index + 2});
        vertices.insert(vertices.end(), temp_vertices.begin(), temp_vertices.end());
    }
    return {vertices, triangles};
}


// TEST(GlobalTest, IGLBunny)
// {
//     init_tables();

//     int nx = 50;
//     int ny = 50;
//     int nz = 50;

//     Eigen::MatrixXd V;
//     Eigen::MatrixXi F;
//     std::string off_path = "D:\\Library\\libigl\\build\\_deps\\libigl_tutorial_tata-src\\bunny.off";
//     igl::readOFF(off_path, V, F);
//     Eigen::Vector3d m = V.colwise().minCoeff();
//     std::cout << m << std::endl;
//     Eigen::Vector3d M = V.colwise().maxCoeff();
//     std::cout << M << std::endl;
//     Eigen::MatrixXd P;
// 	P.resize(nx * ny * nz, 3);
// 	for(unsigned int k = 0; k < nz; k++)
// 	{
// 		double z_axis = m.z() + double(k) * (M.z() - m.z()) / double(nz - 1);
// 		for(unsigned int j = 0; j < ny; j++)
// 		{
// 			double y_axis = m.y() + double(j) * (M.y() - m.y()) / double(ny - 1);
// 			for(unsigned int i = 0; i < nx; i++)
// 			{
// 				double x_axis = m.x() + double(i) * (M.x() - m.x()) / double(nx - 1);
// 				unsigned int count = k * ny * nx + j * nx + i;
// 				P(count, 0) = x_axis;
// 				P(count, 1) = y_axis;
// 				P(count, 2) = z_axis;
// 			}
// 		}
// 	}

//     Eigen::VectorXi I;
//     Eigen::MatrixXd N,C;
//     Eigen::VectorXd S;
//     igl::SignedDistanceType type = igl::SIGNED_DISTANCE_TYPE_PSEUDONORMAL;
//     igl::signed_distance(P, V, F, type, S, I, C, N);

//     // std::cout << "------------------S---------------------" << std::endl;
//     // std::cout << S << std::endl;

//     std::vector<Eigen::Vector3d> vertices;
//     std::vector<Eigen::Vector3i> triangles;
//     int iso_value = 0;

//     for(int k = 0; k < nz - 1; k++)
// 	{
// 		for(int j = 0; j < ny - 1; j++)
// 		{
// 			for(int i = 0; i < nx - 1; i++)
// 			{
// 				//  get cube
//                 std::vector<Eigen::Vector3d> coors{
//                     {static_cast<double>(i), static_cast<double>(j), static_cast<double>(k)},
//                     {static_cast<double>(i), static_cast<double>(j + 1), static_cast<double>(k)},
//                     {static_cast<double>(i), static_cast<double>(j + 1), static_cast<double>(k + 1)},
//                     {static_cast<double>(i), static_cast<double>(j), static_cast<double>(k + 1)},
//                     {static_cast<double>(i + 1), static_cast<double>(j), static_cast<double>(k)},
//                     {static_cast<double>(i + 1), static_cast<double>(j + 1), static_cast<double>(k)},
//                     {static_cast<double>(i + 1), static_cast<double>(j + 1), static_cast<double>(k + 1)},
//                     {static_cast<double>(i + 1), static_cast<double>(j), static_cast<double>(k + 1)},
//                 };
// 				unsigned int count = k * ny * nx + j * nx + i;
//                 std::vector<double> signed_distance{
//                     S[k * ny * nx + j * nx + i],
//                     S[k * ny * nx + (j + 1) * nx + i],
//                     S[(k + 1) * ny * nx + (j + 1) * nx + i],
//                     S[(k + 1) * ny * nx + j * nx + i],
//                     S[k * ny * nx + j * nx + i + 1],
//                     S[k * ny * nx + (j + 1) * nx + i + 1],
//                     S[(k + 1) * ny * nx + (j + 1) * nx + i + 1],
//                     S[(k + 1) * ny * nx + j * nx + i + 1],
//                 };

//                 std::bitset<8> distance_symbol;
//                 for(auto distance = signed_distance.begin(); distance != signed_distance.end(); distance++)
//                 {
//                     int index = std::distance(signed_distance.begin(), distance);
//                     if(*distance - iso_value < 0)
//                     {
//                         distance_symbol.set(index, false);
//                     }
//                     else
//                     {
//                         distance_symbol.set(index, true);
//                     }
//                 }
//                 std::vector<unsigned short int> edgess;
//                 try
//                 {
//                     // edgess = MC_TABLES.at(distance_symbol);
//                     edgess = get_edgess(signed_distance);
//                 }
//                 catch (const std::out_of_range& e)
//                 {
//                     std::cerr << "Caught std::out_of_range exception: " << e.what() << std::endl;
//                     throw e;
//                 }
//                 for(auto edges: edgess)
//                 {
//                     Eigen::Vector3i triangle;
//                     int triangle_edge_count = 0;
//                     while(triangle_edge_count <= 2)
//                     {
//                         auto edge = static_cast<Edge>(edges & 0xF);
//                         auto vertex = get_vertex_by_edge(edge, coors, signed_distance);
//                         vertices.emplace_back(std::move(vertex));
//                         triangle[triangle_edge_count] = vertices.size() - 1;
//                         triangle_edge_count++;
//                         edges = edges >> 4;
//                     }
//                     triangles.emplace_back(std::move(triangle));
//                 }
//             }
//         }
//     }
//     write_obj("./bunny.obj", vertices, triangles);
// }

