#ifndef __FEATURE_POINT_MC33_TABLES_H__
#define __FEATURE_POINT_MC33_TABLES_H__

#include <vector>
#include <unordered_map>
#include <bitset>
#include <unordered_set>
#include <format>
#include <functional>

#include <boost/functional/hash.hpp>
#include <Eigen/Dense>

#include "charles_mc33_type.h"

class FeatureMC33Table
{
public:
    // example: {{4, 0x0387}}, represents feature point interpolated by edge {0, 3, 8, 7}, first element represents length, in case of first is zero
    std::vector<std::tuple<unsigned short, unsigned int>> feature_interpolation_rules;
    // example: {{0, {0x38, 0x37}}}, represents feature point index 0 (0x038) + two point on edge, makes a triangle
    std::unordered_map<unsigned short, std::vector<unsigned short>> fp_connected_edges;
    // with mc33 inserted point triangles, such as{0x38}
    std::vector<unsigned short> mc33_triangles;
    // example: {0x038}, represents normal triangles
    std::vector<unsigned short> common_triangles;
    /**
     * @brief feature point constrains
     * this is mesh(composed by vertex of cube, can be saw as triangles, use with cube coors(can be saw as vertices))
     * 
     * for example:
     * {
     *     {
     *         0,  // this is feature point index
     *         {
     *             {1, 2, 3},  // this is mesh triangle(composed by vertex of cube)
     *             {1, 2, 6},
     *             {2, 3, 6},
     *             {Vertex::vc, 3, 6},  // Vertex::vc is the center point of cube
     *             {Vertex::vc, 1, 6},
     *             {Vertex::vc, 1, 3},
     *         }
     *     }
     * }
     */
    std::unordered_map<unsigned short, std::vector<Eigen::Vector3i>> constrains;
    static int feature_triangle_length;
    static int mc33_triangle_length;
    static int common_triangle_length;

    FeatureMC33Table rotate(const Axis& axis, int times = 1);

    std::unordered_map<unsigned short, std::vector<Eigen::Vector3i>> constrains_rotate(const Axis& axis, int times = 1);

    friend auto operator<<(std::ostream& os, FeatureMC33Table const& m) -> std::ostream&
    {
        // return os << std::format("x = {}, y = {}, z = {}, flag = {}", m.x, m.y, m.z, m.flag);
        os << "feature interpolation rules: {";
        for(const auto& [key, value]: m.feature_interpolation_rules)
        {
            os << std::format(", length: {}, feature point {:#x} ", key, value);
        }
        os << "}, feature triangles: {";
        for(const auto& [key, value]: m.fp_connected_edges)
        {
            os << std::format("feature index: {}, ", key);
            os << "feature edge: {";
            for(const auto& feature_edge: value)
            {
                os << std::format("{:#x}, ", feature_edge);
            }
        }
        os << "}, mc33 triangles: {";
        for(const auto& mc33_triangle: m.mc33_triangles)
        {
            os << std::format("{:#x}, ", mc33_triangle);
        }
        os << "}, common triangles: {";
        for(const auto& common_triangle: m.common_triangles)
        {
            os << std::format("{:#x}, ", common_triangle);
        }
        os << "}, constrains: {";
        for(const auto& [feature_point_index, triangles]: m.constrains)
        {
            os << feature_point_index << "{";
            for(const auto& triangle: triangles)
            {
                os << "{" << triangle[0] << ", " << triangle[1] << ", " << triangle[2] << "}, ";
            }
            os << "}";
        }
        os << "}";
        return os;
    }

    // std::tuple<std::vector<Eigen::Vector3d>, std::vector<Eigen::Vector3i>> get_vertices_triangles(const std::vector<int>& signed_distances, const int& iso_value)
    // {
    //     std::bitset<8> distance_signs;
    //     for(auto signed_distance = signed_distances.begin(); signed_distance != signed_distances.end(); signed_distance++)
    //     {
    //         int index = std::distance(signed_distances.begin(), signed_distance);
    //         if(*signed_distance - iso_value < 0)
    //         {
    //             distance_signs.set(index, false);
    //         }
    //         else
    //         {
    //             distance_signs.set(index, true);
    //         }
    //     }
    // }
};

std::ostream& operator<<(std::ostream& stream, const std::unordered_set<std::unordered_set<Vertex>, boost::hash<std::unordered_set<Vertex>>>& object);

extern std::unordered_map<
    // distance signs
    std::bitset<8>,
    std::unordered_map<
        // connected vertices
        std::unordered_set<std::unordered_set<Vertex>, boost::hash<std::unordered_set<Vertex>>>,
        // triangles(may be this should be reimplemented by class, to make it easier)
        FeatureMC33Table,
        boost::hash<std::unordered_set<std::unordered_set<Vertex>, boost::hash<std::unordered_set<Vertex>>>>
    >
> MC33_TABLES;

extern std::unordered_map<std::unordered_set<Vertex>, Edge, boost::hash<std::unordered_set<Vertex>>> VERTEX_EDGE;

extern std::unordered_map<Edge, std::unordered_set<Vertex>> EDGE_VERTEX;

bool vertex_connected(const std::bitset<8>& distance_sign, const Vertex &vertex1, const Vertex &vertex2);

bool vertex_interpolation_connected(const std::vector<double>& signed_distance, const Vertex& vertex1, const Vertex& vertex2);

void init_tables();


template<typename R, typename... Args>
std::vector<R> handle_edges(const unsigned short& in_edges, const int& edge_num, std::function<R(const Edge&, const Args&...)> callback, const Args&... args)
{
    std::vector<R> ret;
    unsigned short edges = in_edges;
    for(unsigned short i = 0; i < edge_num; i++)
    {
        edges = edges >> (i == 0 ? 0 : 4);
        auto edge = edges & 0xF;
        R res = callback(static_cast<Edge>(edge), args...);
        ret.emplace_back(res);
    }
    return ret;
}

template<typename... Args>
unsigned short handle_edges(const unsigned short& in_edges, const int& edge_num, std::function<Edge(const Edge&, const Args&...)> callback, const Args&... args)
{
    unsigned short res_edges = 0x00000000;
    unsigned short edges = in_edges;
    for(unsigned short i = 0; i < edge_num; i++)
    {
        edges = edges >> (i == 0 ? 0 : 4);
        auto edge = edges & 0xF;
        Edge res_edge = callback(static_cast<Edge>(edge), args...);
        res_edges = (res_edges | (static_cast<unsigned short>(res_edge) << (4 * i)));
    }
    return res_edges;
}

void print_mc33_table();

#endif