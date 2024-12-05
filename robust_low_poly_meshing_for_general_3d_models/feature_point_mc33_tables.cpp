#include "feature_point_mc33_tables.h"
#include "gtest/gtest.h"

#include <unordered_map>
#include <unordered_set>
#include <boost/functional/hash.hpp>
#include <bitset>
#include <queue>
#include <stdexcept>
#include <format>
#include "Eigen/Dense"

#include "charles_mc33_type.h"
#include "tuple_hash.h"

/*
  Vertices:            Edges:               Faces:
    3 ___________2        _____2______         ____________
   /|           /|      /|           /|      /|           /|
  / |          / |     B |          A |     / |    2     / |
7/___________6/  |    /_____6_____ /  |    /___________ /  |
|   |        |   |   |   3        |   1   |   |     4  |   |
|   |        |   |   |   |        |   |   | 3 |        | 1 |     z
|   0________|___1   |   |_____0__|___|   |   |_5______|___|     |
|  /         |  /    7  /         5  /    |  /         |  /      |____y
| /          | /     | 8          | 9     | /      0   | /      /
4/___________5/      |/_____4_____|/      |/___________|/      x

*/

int FeatureMC33Table::feature_triangle_length = 2;
int FeatureMC33Table::mc33_triangle_length = 2;
int FeatureMC33Table::common_triangle_length = 3;

std::unordered_map<Vertex, std::unordered_set<Vertex>> VERTEX_NEIGHBOOR = {
    {Vertex::v0, {Vertex::v1, Vertex::v3, Vertex::v4}},
    {Vertex::v1, {Vertex::v0, Vertex::v2, Vertex::v5}},
    {Vertex::v2, {Vertex::v1, Vertex::v3, Vertex::v6}},
    {Vertex::v3, {Vertex::v0, Vertex::v2, Vertex::v7}},
    {Vertex::v4, {Vertex::v0, Vertex::v5, Vertex::v7}},
    {Vertex::v5, {Vertex::v1, Vertex::v4, Vertex::v6}},
    {Vertex::v6, {Vertex::v2, Vertex::v5, Vertex::v7}},
    {Vertex::v7, {Vertex::v3, Vertex::v4, Vertex::v6}},
};

std::unordered_map<std::unordered_set<Vertex>, Edge, boost::hash<std::unordered_set<Vertex>>> VERTEX_EDGE = {
    {{Vertex::v0, Vertex::v1}, Edge::e0},
    {{Vertex::v1, Vertex::v2}, Edge::e1},
    {{Vertex::v2, Vertex::v3}, Edge::e2},
    {{Vertex::v3, Vertex::v0}, Edge::e3},
    {{Vertex::v4, Vertex::v5}, Edge::e4},
    {{Vertex::v5, Vertex::v6}, Edge::e5},
    {{Vertex::v6, Vertex::v7}, Edge::e6},
    {{Vertex::v7, Vertex::v4}, Edge::e7},
    {{Vertex::v0, Vertex::v4}, Edge::e8},
    {{Vertex::v1, Vertex::v5}, Edge::e9},
    {{Vertex::v2, Vertex::v6}, Edge::e10},
    {{Vertex::v3, Vertex::v7}, Edge::e11},
};

std::unordered_map<Edge, std::unordered_set<Vertex>> EDGE_VERTEX = {
    {Edge::e0, {Vertex::v0, Vertex::v1}},
    {Edge::e1, {Vertex::v1, Vertex::v2}},
    {Edge::e2, {Vertex::v2, Vertex::v3}},
    {Edge::e3, {Vertex::v3, Vertex::v0}},
    {Edge::e4, {Vertex::v4, Vertex::v5}},
    {Edge::e5, {Vertex::v5, Vertex::v6}},
    {Edge::e6, {Vertex::v6, Vertex::v7}},
    {Edge::e7, {Vertex::v7, Vertex::v4}},
    {Edge::e8, {Vertex::v0, Vertex::v4}},
    {Edge::e9, {Vertex::v1, Vertex::v5}},
    {Edge::e10, {Vertex::v2, Vertex::v6}},
    {Edge::e11, {Vertex::v3, Vertex::v7}},
};

std::unordered_map<std::unordered_set<Vertex>, unsigned short, boost::hash<std::unordered_set<Vertex>>> NEIGHBOR_DISTANCE = {
    {{Vertex::v0, Vertex::v1}, 1},
    {{Vertex::v0, Vertex::v2}, 2},
    {{Vertex::v0, Vertex::v3}, 1},
    {{Vertex::v0, Vertex::v4}, 1},
    {{Vertex::v0, Vertex::v5}, 2},
    {{Vertex::v0, Vertex::v6}, 3},
    {{Vertex::v0, Vertex::v7}, 2},

    {{Vertex::v1, Vertex::v2}, 1},
    {{Vertex::v1, Vertex::v3}, 2},
    {{Vertex::v1, Vertex::v4}, 2},
    {{Vertex::v1, Vertex::v5}, 1},
    {{Vertex::v1, Vertex::v6}, 2},
    {{Vertex::v1, Vertex::v7}, 3},

    {{Vertex::v2, Vertex::v3}, 1},
    {{Vertex::v2, Vertex::v4}, 3},
    {{Vertex::v2, Vertex::v5}, 2},
    {{Vertex::v2, Vertex::v6}, 1},
    {{Vertex::v2, Vertex::v7}, 2},

    {{Vertex::v3, Vertex::v4}, 2},
    {{Vertex::v3, Vertex::v5}, 3},
    {{Vertex::v3, Vertex::v6}, 2},
    {{Vertex::v3, Vertex::v7}, 1},

    {{Vertex::v4, Vertex::v5}, 1},
    {{Vertex::v4, Vertex::v6}, 2},
    {{Vertex::v4, Vertex::v7}, 1},

    {{Vertex::v5, Vertex::v6}, 1},
    {{Vertex::v5, Vertex::v7}, 2},

    {{Vertex::v6, Vertex::v7}, 1},
};

std::unordered_map<std::unordered_set<Vertex>, Face, boost::hash<std::unordered_set<Vertex>>> DIAGNAL_VERTEX_FACE = {
    {{Vertex::v0, Vertex::v5}, Face::f0},
    {{Vertex::v1, Vertex::v4}, Face::f0},
    {{Vertex::v1, Vertex::v6}, Face::f1},
    {{Vertex::v2, Vertex::v5}, Face::f1},
    {{Vertex::v2, Vertex::v7}, Face::f2},
    {{Vertex::v3, Vertex::v6}, Face::f2},
    {{Vertex::v0, Vertex::v7}, Face::f3},
    {{Vertex::v3, Vertex::v4}, Face::f3},
    {{Vertex::v0, Vertex::v2}, Face::f4},
    {{Vertex::v1, Vertex::v3}, Face::f4},
    {{Vertex::v4, Vertex::v6}, Face::f5},
    {{Vertex::v5, Vertex::v7}, Face::f5},
};

std::unordered_map<Face, std::unordered_set<Vertex>> FACE_VERTICES = {
    {Face::f0, {Vertex::v0, Vertex::v1, Vertex::v5, Vertex::v4}},
    {Face::f1, {Vertex::v1, Vertex::v2, Vertex::v6, Vertex::v5}},
    {Face::f2, {Vertex::v2, Vertex::v3, Vertex::v7, Vertex::v6}},
    {Face::f3, {Vertex::v0, Vertex::v3, Vertex::v7, Vertex::v4}},
    {Face::f4, {Vertex::v0, Vertex::v1, Vertex::v2, Vertex::v3}},
    {Face::f5, {Vertex::v4, Vertex::v5, Vertex::v6, Vertex::v7}},
};

Face TOP_FACE = Face::f2;

Face BOTTOM_FACE = Face::f0;

std::unordered_map<std::unordered_set<Vertex>, std::unordered_map<std::string, Vertex>, boost::hash<std::unordered_set<Vertex>>> CUBE_INTERPOLATION_MAP = {
    {
        {Vertex::v0, Vertex::v6},
        {
            {"A0", Vertex::v0},
            {"C0", Vertex::v5},
            {"B0", Vertex::v1},
            {"D0", Vertex::v4},
            {"A1", Vertex::v3},
            {"C1", Vertex::v6},
            {"B1", Vertex::v2},
            {"D1", Vertex::v7}
        }
    },
    {
        {Vertex::v1, Vertex::v7},
        {
            {"A0", Vertex::v1},
            {"C0", Vertex::v4},
            {"B0", Vertex::v5},
            {"D0", Vertex::v0},
            {"A1", Vertex::v2},
            {"C1", Vertex::v7},
            {"B1", Vertex::v6},
            {"D1", Vertex::v3}
        }
    },
    {
        {Vertex::v5, Vertex::v3},
        {
            {"A0", Vertex::v5},
            {"C0", Vertex::v0},
            {"B0", Vertex::v4},
            {"D0", Vertex::v1},
            {"A1", Vertex::v6},
            {"C1", Vertex::v3},
            {"B1", Vertex::v7},
            {"D1", Vertex::v2}
        }
    },
    {
        {Vertex::v4, Vertex::v2},
        {
            {"A0", Vertex::v4},
            {"C0", Vertex::v1},
            {"B0", Vertex::v0},
            {"D0", Vertex::v5},
            {"A1", Vertex::v7},
            {"C1", Vertex::v2},
            {"B1", Vertex::v3},
            {"D1", Vertex::v6}
        }
    }
};

// clock wise
std::unordered_map<Vertex, Vertex> VERTEX_X_ROTATE = {
    {Vertex::v0, Vertex::v3},
    {Vertex::v1, Vertex::v0},
    {Vertex::v2, Vertex::v1},
    {Vertex::v3, Vertex::v2},
    {Vertex::v4, Vertex::v7},
    {Vertex::v5, Vertex::v4},
    {Vertex::v6, Vertex::v5},
    {Vertex::v7, Vertex::v6},
    {Vertex::vc, Vertex::vc},
};

// clock wise
std::unordered_map<Vertex, Vertex> VERTEX_Y_ROTATE = {
    {Vertex::v0, Vertex::v4},
    {Vertex::v1, Vertex::v5},
    {Vertex::v2, Vertex::v1},
    {Vertex::v3, Vertex::v0},
    {Vertex::v4, Vertex::v7},
    {Vertex::v5, Vertex::v6},
    {Vertex::v6, Vertex::v2},
    {Vertex::v7, Vertex::v3},
    {Vertex::vc, Vertex::vc},
};

// clock wise
std::unordered_map<Vertex, Vertex> VERTEX_Z_ROTATE = {
    {Vertex::v0, Vertex::v1},
    {Vertex::v1, Vertex::v5},
    {Vertex::v2, Vertex::v6},
    {Vertex::v3, Vertex::v2},
    {Vertex::v4, Vertex::v0},
    {Vertex::v5, Vertex::v4},
    {Vertex::v6, Vertex::v7},
    {Vertex::v7, Vertex::v3},
    {Vertex::vc, Vertex::vc},
};


bool vertex_in_top_face(const Vertex& vertex)
{
    auto vertices = FACE_VERTICES.at(TOP_FACE);
    if(vertices.contains(vertex))
    {
        return true;
    }
    return false;
}

bool vertex_in_bottom_face(const Vertex& vertex)
{
    auto vertices = FACE_VERTICES.at(BOTTOM_FACE);
    if(vertices.contains(vertex))
    {
        return true;
    }
    return false;
}

bool has_value_bigger_than_zero_in_interval(
    const double& a, const double& b, const double& c,
    const double& start = 0.0f, const double& end = 1.0f
)
{
    auto fx = [](double a, double b, double c, double x) -> double{
        return a * std::pow(x, 2) + b * x + c;
    };
    if(a >= 0.0f)
    {
        if(fx(a, b, c, start) > 0 || fx(a, b, c, end) > 0)
        {
            return true;
        }
        return false;
    }
    else
    {
        // must has root, so(b^2 - 4ac > 0)
        if(std::pow(b, 2) - 4 * a * c > 0)
        {
            // if (fx(start) > 0) or (fx(end) > 0) or (-b/2a between start and end)
            if(fx(a, b, c, start) > 0 || fx(a, b, c, end) > 0 || (-b / (2.0f * a) > start && -b / (2.0f * a) < end))
            {
                return true;
            }
            return false;
        }
        return false;
    }
}

/**
 * @brief suppose vertex1 and vertex2 are not obviously connected
 * 
 * @param signed_distance 
 * @param vertex1 
 * @param vertex2 
 * @return true 
 * @return false 
 */
bool vertex_interpolation_connected(const std::vector<double>& signed_distance, const Vertex& vertex1, const Vertex& vertex2)
{
    // only has two cases, one is in the diagonal in same face, one is in diagonal of the cube
    // check if two vertex is one the same face(easy, check distance will be fine)
    auto distance = NEIGHBOR_DISTANCE.at(std::unordered_set<Vertex>{vertex1, vertex2});
    switch(distance)
    {
        case 1:
            // distance is 1, is same edge, must be connected
            return true;
        case 2:
        {
            // distance is 2, is diagnal point in same face, consider interpolation connected
            auto face = DIAGNAL_VERTEX_FACE.at(std::unordered_set<Vertex>{vertex1, vertex2});
            auto other_two_points = FACE_VERTICES.at(face);
            other_two_points.erase(vertex1);
            other_two_points.erase(vertex2);
            double distance_ac = signed_distance[static_cast<unsigned int>(vertex1)] * signed_distance[static_cast<unsigned int>(vertex2)];
            double distance_bd = 1.0f;
            for(const auto& elem: other_two_points)
            {
                distance_bd = distance_bd * signed_distance[static_cast<unsigned int>(elem)];
            }
            if(distance_ac - distance_bd > 0.0f)
            {
                return true;
            }
            break;
        }
        case 3:
        {
            // distance is 3, is diagnal point in cube, just compute the top and bottom interpolation
            auto cube_interpolation = CUBE_INTERPOLATION_MAP.at(std::unordered_set<Vertex>{vertex1, vertex2});
            double distance_a0 = signed_distance[static_cast<unsigned short>(cube_interpolation["A0"])];
            double distance_a1 = signed_distance[static_cast<unsigned short>(cube_interpolation["A1"])];
            double distance_b0 = signed_distance[static_cast<unsigned short>(cube_interpolation["B0"])];
            double distance_b1 = signed_distance[static_cast<unsigned short>(cube_interpolation["B1"])];
            double distance_c0 = signed_distance[static_cast<unsigned short>(cube_interpolation["C0"])];
            double distance_c1 = signed_distance[static_cast<unsigned short>(cube_interpolation["C1"])];
            double distance_d0 = signed_distance[static_cast<unsigned short>(cube_interpolation["D0"])];
            double distance_d1 = signed_distance[static_cast<unsigned short>(cube_interpolation["D1"])];
            auto a = (distance_a1 - distance_a0)*(distance_c1 - distance_c0) - (distance_b1 - distance_b0) * (distance_d1 - distance_d0);
            auto b = distance_c0*(distance_a1 - distance_a0) + distance_a0*(distance_c1 - distance_c0) - distance_d0*(distance_b1 - distance_b0) - distance_b0*(distance_d1 - distance_d0);
            auto c = distance_a0*distance_c0 - distance_b0 * distance_d0;
            if(has_value_bigger_than_zero_in_interval(a, b, c, 0.0f, 1.0f))
            {
                return true;
            }
            // if(a < 0)
            // {
            //     auto t_max = -b / (2 * a);
            //     if(t_max > 0.0f && t_max < 1.0f)
            //     {
            //         if(a * std::pow(t_max, 2) + b * t_max + c > 0)
            //         {
            //             // joined
            //             return true;
            //         }
            //     }
            // }
            break;
        }
        default:
            break;
    }
    return false;
}

std::unordered_map<std::tuple<std::bitset<8>, Vertex, Vertex>, bool, hash_tuple> CONNECTED_VERTEX_CACHE;

/**
 * @brief find if there are path between v1 and v2
 * 
 * @param vertex1 
 * @param vertex2 
 * @return true 
 * @return false 
 */
bool vertex_connected(const std::bitset<8>& distance_sign, const Vertex &vertex1, const Vertex &vertex2)
{
    try
    {
        auto is_connected = CONNECTED_VERTEX_CACHE.at({distance_sign, vertex1, vertex2});
        return is_connected;
    }
    catch(const std::out_of_range& e)
    {
        // do nothing
    }
    auto set_cache_lambda = [&](bool is_connected){
        CONNECTED_VERTEX_CACHE.emplace(std::make_tuple(distance_sign, vertex1, vertex2), is_connected);
    };
    if(vertex1 == vertex2)
    {
        set_cache_lambda(true);
        return true;
    }
    auto vertex1_symbol = distance_sign[static_cast<unsigned short>(vertex1)];
    auto vertex2_symbol = distance_sign[static_cast<unsigned short>(vertex2)];
    if(vertex1_symbol != vertex2_symbol)
    {
        set_cache_lambda(false);
        return false;
    }
    auto vertex_symbol = vertex1_symbol;

    // avoid visited vertex
    std::unordered_set<Vertex> visited_vertices;
    std::queue<std::pair<Vertex, unsigned short>> visiting_vertices;
    visiting_vertices.push(std::make_pair(vertex1, 0));
    while (!visiting_vertices.empty())
    {
        auto current_vertex = visiting_vertices.front();
        visiting_vertices.pop();
        auto depth = current_vertex.second;
        if (depth >= 4)
        {
            set_cache_lambda(false);
            return false;
        }
        if (current_vertex.first == vertex2)
        {
            set_cache_lambda(true);
            return true;
        }
        auto neighboors = VERTEX_NEIGHBOOR.at(current_vertex.first);
        for(const auto& neighboor: neighboors)
        {
            // symbol is same
            if(distance_sign[static_cast<unsigned short>(neighboor)] == vertex_symbol)
            {
                if(!visited_vertices.contains(neighboor))
                {
                    visiting_vertices.push(std::make_pair(neighboor, depth+1));
                }
            }
        }
        visited_vertices.emplace(current_vertex.first);
    }
    set_cache_lambda(false);
    return false;
}

Vertex rotate_vertex_in_axis(const Vertex& vertex, const Axis& axis, int times = 1)
{
    Vertex rotated_vertex = vertex;
    for(int i = 0; i < times; i++)
    {
        switch(axis)
        {
        case Axis::x:
            rotated_vertex = VERTEX_X_ROTATE.at(rotated_vertex);
            break;
        case Axis::y:
            rotated_vertex = VERTEX_Y_ROTATE.at(rotated_vertex);
            break;
        case Axis::z:
            rotated_vertex = VERTEX_Z_ROTATE.at(rotated_vertex);
            break;
        default:
            throw std::invalid_argument("axis is weird");
            break;
        }
    }
    return rotated_vertex;
}

std::vector<Vertex> rotate_vertices_in_axis(const std::vector<Vertex>& vertices, const Axis& axis, int times = 1)
{
    std::vector<Vertex> rotated_vertices;
    for(const auto& vertex: vertices)
    {
        auto rotated_vertex = vertex;
        rotated_vertex = rotate_vertex_in_axis(rotated_vertex, axis, times);
        rotated_vertices.emplace_back(rotated_vertex);
    }
    return rotated_vertices;
}

Edge rotate_edge_in_axis(const Edge& edge, const Axis& axis, int times = 1)
{
    auto vertices = EDGE_VERTEX.at(edge);
    std::unordered_set<Vertex> rotated_vertices;
    for(const auto& vertex: vertices)
    {
        auto rotated_vertex = vertex;
        for(int i = 0; i < times; i++)
        {
            switch (axis)
            {
            case Axis::x:
                rotated_vertex = VERTEX_X_ROTATE.at(rotated_vertex);
                break;
            case Axis::y:
                rotated_vertex = VERTEX_Y_ROTATE.at(rotated_vertex);
                break;
            case Axis::z:
                rotated_vertex = VERTEX_Z_ROTATE.at(rotated_vertex);
                break;
            default:
                throw std::invalid_argument("axis is weird");
                break;
            }
        }
        rotated_vertices.emplace(rotated_vertex);
    }
    auto rotated_edge = VERTEX_EDGE.at(rotated_vertices);
    return rotated_edge;
}

unsigned short handle_edges(const unsigned short& in_edges, const int& edge_num, std::function<Edge(const Edge&)> callback)
{
    unsigned short res_edges = 0x00000000;
    unsigned short edges = in_edges;
    for(unsigned short i = 0; i < edge_num; i++)
    {
        edges = edges >> (i == 0 ? 0 : 4);
        auto edge = edges & 0xF;
        Edge res_edge = callback(static_cast<Edge>(edge));
        res_edges = (res_edges | (static_cast<unsigned short>(res_edge) << (4 * i)));
    }
    return res_edges;
}

unsigned short rotate_edge_in_axis(unsigned short edges, const Axis& axis, int times = 1, int edge_num = 3)
{
    unsigned short rotated_edges = 0x00000000;
    for(unsigned short i = 0; i < edge_num; i++)
    {
        edges = edges >> (i == 0 ? 0 : 4);
        auto edge = edges & 0xF;
        Edge rotated_edge;
        rotated_edge = rotate_edge_in_axis(static_cast<Edge>(edge), axis, times);
        rotated_edges = (rotated_edges | (static_cast<unsigned short>(rotated_edge) << (4 * i)));
    }
    return rotated_edges;
}

std::vector<unsigned short> rotate_edge_in_axis(const std::vector<unsigned short>& edgess, const Axis& axis, int times = 1, int length = 3)
{
    std::vector<unsigned short> rotated_edgess;
    for(auto edges: edgess)
    {
        unsigned short rotated_edges = rotate_edge_in_axis(edges, axis, times, length);
        // unsigned short rotated_edges = 0x0000;
        // for(int i = 0; i < 3; i++)
        // {
        //     edges = edges >> (i == 0 ? 0 : 4);
        //     auto edge = edges & 0xF;
        //     Edge rotated_edge;
        //     rotated_edge = rotate_edge_in_axis(static_cast<Edge>(edge), axis, times);
        //     rotated_edges = (rotated_edges | (static_cast<unsigned short>(rotated_edge) << (4 * i)));
        // }
        rotated_edgess.emplace_back(rotated_edges);
    }
    return rotated_edgess;
}

std::vector<Vertex> get_negative_vertices_from_signed_distance(const std::bitset<8>& signed_distance)
{
    std::vector<Vertex> vertices;
    for(int i = 0; i < signed_distance.size(); i++)
    {
        if(signed_distance[i] == false)
        {
            vertices.emplace_back(static_cast<Vertex>(i));
        }
    }
    return vertices;
}

std::unordered_set<std::unordered_set<Vertex>, boost::hash<std::unordered_set<Vertex>>> rotate_vertices_in_axis(
    const std::unordered_set<std::unordered_set<Vertex>, boost::hash<std::unordered_set<Vertex>>>& verticess,
    const Axis& axis, int times = 1
)
{
    std::unordered_set<std::unordered_set<Vertex>, boost::hash<std::unordered_set<Vertex>>> rotated_verticess;
    for(const auto& connected_vertices: verticess)
    {
        std::unordered_set<Vertex> rotated_connnected_vertices;
        for(const auto& connected_vertex: connected_vertices)
        {
            Vertex rotated_connected_vertex = rotate_vertex_in_axis(connected_vertex, axis, times);
            rotated_connnected_vertices.emplace(rotated_connected_vertex);
        }
        rotated_verticess.emplace(rotated_connnected_vertices);
    }
    return rotated_verticess;
}

void init_tables()
{
    // rotate table
    std::unordered_map<
        std::bitset<8>,
        std::unordered_map<
            std::unordered_set<std::unordered_set<Vertex>, boost::hash<std::unordered_set<Vertex>>>,
            FeatureMC33Table,
            boost::hash<std::unordered_set<std::unordered_set<Vertex>, boost::hash<std::unordered_set<Vertex>>>>
        >
    > rotated_tables;
    for(auto [distance_signs, sub_cases]: MC33_TABLES)
    {
        if (distance_signs == 0b11011100)
        {
            std::cout << "catch point" << std::endl;
        }
        // get negative vertex
        auto vertices = get_negative_vertices_from_signed_distance(distance_signs);
        for(auto [connected_verticess, feature_mc33_table]: sub_cases)
        {
            for(unsigned short k = 0; k < 4; k++)
            {
                auto z_rotated_vertices = rotate_vertices_in_axis(vertices, Axis::z, k);
                // auto z_rotated_edgess = rotate_edge_in_axis(triangles, Axis::z, k);
                auto z_rotated_table = feature_mc33_table.rotate(Axis::z, k);
                auto z_rotated_connected_vertices = rotate_vertices_in_axis(connected_verticess, Axis::z, k);
                for(unsigned short j = 0; j < 4; j++)
                {
                    auto y_rotated_vertices = rotate_vertices_in_axis(z_rotated_vertices, Axis::y, j);
                    // auto y_rotated_edgess = rotate_edge_in_axis(z_rotated_edgess, Axis::y, j);
                    auto y_rotated_table = z_rotated_table.rotate(Axis::y, j);
                    auto y_rotated_connected_vertices = rotate_vertices_in_axis(z_rotated_connected_vertices, Axis::y, j);
                    for(unsigned short i = 0; i < 4; i++)
                    {
                        auto x_rotated_vertices = rotate_vertices_in_axis(y_rotated_vertices, Axis::x, i);
                        // auto x_rotated_edgess = rotate_edge_in_axis(y_rotated_edgess, Axis::x, i);
                        auto x_rotated_table = y_rotated_table.rotate(Axis::x, i);
                        auto x_rotated_connected_vertices = rotate_vertices_in_axis(y_rotated_connected_vertices, Axis::x, i);
                        // use rotated vertices to construct new signed distance
                        std::bitset<8> rotated_signed_distance{0b11111111};
                        for(const auto& x_rotated_vertice: x_rotated_vertices)
                        {
                            rotated_signed_distance.set(static_cast<int>(x_rotated_vertice), false);
                        }
                        if(rotated_signed_distance == 0b10011011)
                        {
                            std::cout << "catch point" << std::endl;
                        }
                        if(!rotated_tables.contains(rotated_signed_distance))
                        {
                            std::unordered_map<
                                std::unordered_set<std::unordered_set<Vertex>, boost::hash<std::unordered_set<Vertex>>>,
                                FeatureMC33Table,
                                boost::hash<std::unordered_set<std::unordered_set<Vertex>, boost::hash<std::unordered_set<Vertex>>>>
                            > temp_map;
                            temp_map.emplace(std::make_pair(x_rotated_connected_vertices, x_rotated_table));
                            rotated_tables.emplace(std::make_pair(rotated_signed_distance, temp_map));
                        }
                        else
                        {
                            if(!rotated_tables.at(rotated_signed_distance).contains(x_rotated_connected_vertices))
                            {
                                rotated_tables.at(rotated_signed_distance).emplace(std::make_pair(x_rotated_connected_vertices, x_rotated_table));
                            }
                        }
                    }
                }
            }
        }
    }

    MC33_TABLES.insert(rotated_tables.begin(), rotated_tables.end());
    std::cout << "after rotate, mc33 table size is: " << MC33_TABLES.size() << std::endl;

    // invert distance signs
    std::unordered_map<
        std::bitset<8>,
        std::unordered_map<
            std::unordered_set<std::unordered_set<Vertex>, boost::hash<std::unordered_set<Vertex>>>,
            FeatureMC33Table,
            boost::hash<std::unordered_set<std::unordered_set<Vertex>, boost::hash<std::unordered_set<Vertex>>>>
        >
    > invert_table;
    for(auto [key, value]: MC33_TABLES)
    {
        auto origin = key;
        origin.flip();
        invert_table.emplace(std::make_pair(origin, value));
    }

    MC33_TABLES.insert(invert_table.begin(), invert_table.end());

    std::cout << "after invert, mc33 table size is: " << MC33_TABLES.size() << std::endl;
}

void print_mc33_table()
{
    for(const auto& [distances_signs, sub_cases]: MC33_TABLES)
    {
        std::cout << "{ distances_signs: "<< distances_signs << ", ";
        std::cout << "sub_cases: {";
        for(const auto& [connected_verticess, feature_mc33_table]: sub_cases)
        {
            std::cout << "{" << "connected_vertices: {";
            for(const auto& connected_vertices: connected_verticess)
            {
                std::cout << "{";
                for(const auto& connected_vertex: connected_vertices)
                {
                    std::cout << static_cast<unsigned short>(connected_vertex) << " ";
                }
                std::cout << "},";
            }
            std::cout << "}, mc33 feature table:{";
            std::cout << feature_mc33_table;
            // for(const auto& triangle: triangles)
            // {
            //     std::cout << std::format("{:#x} ", triangle) << ", ";
            // }
            std::cout << "}";
        }
        std::cout << "}," << std::endl;
    }
}

std::unordered_map<unsigned short, std::vector<Eigen::Vector3i>> FeatureMC33Table::constrains_rotate(const Axis& axis, int times)
{
    std::unordered_map<unsigned short, std::vector<Eigen::Vector3i>> ret;
    for(const auto[feature_point_index, triangles]: this->constrains)
    {
        std::vector<Eigen::Vector3i> rotated_triangles;
        for(const auto& triangle: triangles)
        {
            Eigen::Vector3i rotated_triangle;
            rotated_triangle[0] = static_cast<int>(::rotate_vertex_in_axis(static_cast<Vertex>(triangle[0]), axis, times));
            rotated_triangle[1] = static_cast<int>(::rotate_vertex_in_axis(static_cast<Vertex>(triangle[1]), axis, times));
            rotated_triangle[2] = static_cast<int>(::rotate_vertex_in_axis(static_cast<Vertex>(triangle[2]), axis, times));
            rotated_triangles.emplace_back(std::move(rotated_triangle));
        }
        ret.emplace(std::make_pair(feature_point_index, std::move(rotated_triangles)));
    }
    return ret;
}

FeatureMC33Table FeatureMC33Table::rotate(const Axis& axis, int times)
{
    FeatureMC33Table feature_mc33_tables;
    for(const auto& feature_point: this->feature_interpolation_rules)
    {
        auto [length, fp] = feature_point;
        unsigned short new_fp = ::rotate_edge_in_axis(fp, axis, times, length);
        feature_mc33_tables.feature_interpolation_rules.emplace_back(std::make_tuple(length, new_fp));
    }
    for(auto [index, triangles]: this->fp_connected_edges)
    {
        auto new_triangles = ::rotate_edge_in_axis(triangles, axis, times);
        feature_mc33_tables.fp_connected_edges.emplace(index, new_triangles);
    }
    auto new_mc33_triangles = ::rotate_edge_in_axis(this->mc33_triangles, axis, times, 2);
    //feature_mc33_tables.mc33_triangles.emplace_back(new_mc33_triangles);
    feature_mc33_tables.mc33_triangles.insert(feature_mc33_tables.mc33_triangles.end(), new_mc33_triangles.begin(), new_mc33_triangles.end());
    auto new_common_triangles = ::rotate_edge_in_axis(this->common_triangles, axis, times, 3);
    feature_mc33_tables.common_triangles = new_common_triangles;
    feature_mc33_tables.constrains = this->constrains_rotate(axis, times);
    return feature_mc33_tables;
}

const std::unordered_map<Face, std::vector<Eigen::Vector3i>> FACE_TRIANGLES = {
    {
        Face::f0,
        {
            {0, 1, 5},
            {0, 4, 5},
        }
    },
    {
        Face::f1,
        {
            {1, 2, 6},
            {1, 5, 6},
        }
    },
    {
        Face::f2,
        {
            {2, 3, 6},
            {3, 6, 7}
        }
    },
    {
        Face::f3,
        {
            {0, 3, 7},
            {0, 7, 4}
        }
    },
    {
        Face::f4,
        {
            {1, 2, 3},
            {0, 1, 3}
        }
    },
    {
        Face::f5,
        {
            {5, 6, 7},
            {4, 5, 7}
        }
    }
};

const std::vector<Eigen::Vector3i> CUBE_TRIANGLES = {
    FACE_TRIANGLES.at(Face::f0)[0],
    FACE_TRIANGLES.at(Face::f0)[1],
    FACE_TRIANGLES.at(Face::f1)[0],
    FACE_TRIANGLES.at(Face::f1)[1],
    FACE_TRIANGLES.at(Face::f2)[0],
    FACE_TRIANGLES.at(Face::f2)[1],
    FACE_TRIANGLES.at(Face::f3)[0],
    FACE_TRIANGLES.at(Face::f3)[1],
    FACE_TRIANGLES.at(Face::f4)[0],
    FACE_TRIANGLES.at(Face::f4)[1],
    FACE_TRIANGLES.at(Face::f5)[0],
    FACE_TRIANGLES.at(Face::f5)[1],
};

/**
 * @brief first level deside big case, second level deside divided case
 * first level key is distance_sign, second level is connected vertices of ambiguous
 * 
 */
std::unordered_map<
    // distance signs
    std::bitset<8>,
    std::unordered_map<
        // connected vertices
        std::unordered_set<std::unordered_set<Vertex>, boost::hash<std::unordered_set<Vertex>>>,
        // triangles(may be this should be reimplemented by class, to make it easier)
        FeatureMC33Table,
        boost::hash<std::unordered_set<std::unordered_set<Vertex>, boost::hash<std::unordered_set<Vertex>>>>
    >
> MC33_TABLES = { // NOTE: there are problem in vscode, but no problem in visual studio, we quite sure it's problem of vscode
    // case 0
    {
        // distance signs
        0b11111111,
        {
            {
                // connected vertices(unordered_set<std::unordered_set<Vertex>>)
                {},
                // FeatureMC33Table
                {
                    {}, {}, {}, {}, {}
                }
            }
        }
    },

    // case 1
    {
        // distance signs
        0b11111110,
        // cube
        {
            {
                // connected vertices(unordered_set<std::unordered_set<Vertex>>)
                {},
                // FeatureMC33Table
                {
                    // feature interpolation rules
                    {{3, 0x038}},
                    // feature triangles
                    {
                        {0, {0x38, 0x03, 0x80}}
                    },
                    // mc33 triangles
                    {},
                    // common triangles
                    {},
                    // feature point constrains
                    {
                        {
                            0,
                            CUBE_TRIANGLES
                        }
                    },
                }
            }
        }
    },

    // case 2
    {
        // distance signs
        0b11111100,
        // cube
        {
            {
                // connected vertices(unordered_set<std::unordered_set<Vertex>>)
                {},
                // FeatureMC33Table
                {
                    // feature interpolation rules(std::vector<std::tuple<unsigned short, unsigned short>>)
                    {{4, 0x3891}},
                    // feature triangles(std::unordered_map<unsigned short, std::vector<unsigned short>>)
                    {
                        {0, {0x38, 0x89, 0x91, 0x13}}
                    },
                    // mc33 triangles(std::vector<unsigned short>)
                    {},
                    // common triangles(std::vector<unsigned short>)
                    {},
                    // feature point constrains
                    {
                        {
                            0,
                            CUBE_TRIANGLES
                        }
                    },
                }
            }
        }
    },

    // case 3
    {
        // distance signs
        0b10101111,
        // cube
        {
            // case 3.1
            {
                // connected vertices(unordered_set<std::unordered_set<Vertex>>)
                {{Vertex::v5, Vertex::v7}},
                // FeatureMC33Table
                {
                    // feature interpolation rules
                    {{3, 0x478}, {3, 0x56A}},
                    // feature triangles
                    {
                        {0, {0x74, 0x87, 0x48}},
                        {1, {0x56, 0x6A, 0xA5}}
                    },
                    // mc33 triangles
                    {},
                    // common triangles
                    {},
                    // feature point constrains(std::unordered_map<unsigned short, std::vector<Eigen::Vector3i>>)
                    {
                        {
                            0,
                            {
                                {4, 5, 7},
                                {0, 1, 3},
                                {0, 4, 1},
                                {1, 4, 5},
                                {0, 3, 7},
                                {0, 4, 7},
                                {1, 3, 7},
                                {1, 7, 5},
                            }
                        },
                        {
                            1,
                            {
                                {5, 6, 7},
                                {1, 2, 3},
                                {2, 3, 6},
                                {3, 6, 7},
                                {1, 2, 6},
                                {1, 6, 5},
                                {1, 3, 7},
                                {1, 7, 5},
                            }
                        }
                    },
                }
            },
            // case 3.2
            {
                // connected vertices(unordered_set<std::unordered_set<Vertex>>)
                {{Vertex::v4, Vertex::v6}},
                // FeatureMC33Table
                {
                    // feature interpolation rules
                    {{6, 0x87456A}},
                    // feature triangles
                    {
                        {0, {0x87, 0x48, 0x54, 0xA5, 0x6A, 0x76}}
                    },
                    // mc33 triangles
                    {},
                    // common triangles
                    {},
                    // feature point constrains(std::unordered_map<unsigned short, std::vector<Eigen::Vector3i>>)
                    {
                        {
                            0,
                            CUBE_TRIANGLES
                        }
                    },
                }
            }
        }
    },

    // case 4
    {
        // distance signs
        0b11101011,
        // cube
        {
            // case 4.1.1
            {
                // connected vertices(unordered_set<std::unordered_set<Vertex>>)
                {},
                // FeatureMC33Table
                {
                    // feature interpolation rules(std::vector<std::tuple<unsigned short, unsigned short>>)
                    {{3, 0x12A}, {3, 0x478}},
                    // feature triangles(std::unordered_map<unsigned short, std::vector<unsigned short>>)
                    {
                        {0, {0x21, 0xA2, 0x1A}},
                        {1, {0x48, 0x87, 0x74}}
                    },
                    // mc33 triangles(std::vector<unsigned short>)
                    {},
                    // common triangles(std::vector<unsigned short>)
                    {},
                    // feature point constrains(std::unordered_map<unsigned short, std::vector<Eigen::Vector3i>>)
                    {
                        {
                            1,
                            {
                                {4, 5, 7},
                                {0, 1, 3},
                                {0, 4, 1},
                                {1, 4, 5},
                                {0, 3, 7},
                                {0, 4, 7},
                                {1, 3, 7},
                                {1, 7, 5},
                            }
                        },
                        {
                            0,
                            {
                                {5, 6, 7},
                                {1, 2, 3},
                                {2, 3, 6},
                                {3, 6, 7},
                                {1, 2, 6},
                                {1, 6, 5},
                                {1, 3, 7},
                                {1, 7, 5},
                            }
                        }
                    },
                }
            },
            // case 4.1.2
            {
                // connected vertices(unordered_set<std::unordered_set<Vertex>>)
                {{Vertex::v2, Vertex::v4}},
                // FeatureMC33Table
                {
                    // feature interpolation rules(std::vector<std::tuple<unsigned short, unsigned short>>)
                    {},
                    // feature triangles(std::unordered_map<unsigned short, std::vector<unsigned short>>)
                    {},
                    // mc33 triangles(std::vector<unsigned short>)
                    {},
                    // common triangles(std::vector<unsigned short>)
                    {0x27A, 0xA74, 0xA41, 0x148, 0x218, 0x287},
                    // feature point constrains(std::unordered_map<unsigned short, std::vector<Eigen::Vector3i>>)
                    {},
                }
            },
        }
    },

    // case 5
    {
        // distance signs
        0b11011100,
        // cube
        {
            {
                // connected vertices(unordered_set<std::unordered_set<Vertex>>)
                {},
                // FeatureMC33Table
                {
                    // feature interpolation rules(std::vector<std::tuple<unsigned short, unsigned short>>)
                    {{5, 0x13845}},
                    // feature triangles(std::unordered_map<unsigned short, std::vector<unsigned short>>)
                    {
                        {0, {0x13, 0x38, 0x84, 0x45, 0x51}}
                    },
                    // mc33 triangles(std::vector<unsigned short>)
                    {},
                    // common triangles(std::vector<unsigned short>)
                    {},
                    // feature point constrains(std::unordered_map<unsigned short, std::vector<Eigen::Vector3i>>)
                    {
                        {
                            0,
                            CUBE_TRIANGLES
                        }
                    },
                }
            }
        }
    },

    // case 6
    {
        // distance signs
        0b11001011,
        // cube
        {
            // case 6.1.1
            {
                // connected vertices(unordered_set<std::unordered_set<Vertex>>)
                {{Vertex::v1, Vertex::v6}},
                // FeatureMC33Table
                {
                    // feature interpolation rules(std::vector<std::tuple<unsigned short, unsigned short>>)
                    {{4, 0x8957}, {3, 0x12A}},
                    // feature triangles(std::unordered_map<unsigned short, std::vector<unsigned short>>)
                    {
                        {0, {0x98, 0x59, 0x75, 0x87}},
                        {1, {0x21, 0xA2, 0x1A}}
                    },
                    // mc33 triangles(std::vector<unsigned short>)
                    {},
                    // common triangles(std::vector<unsigned short>)
                    {},
                    // feature point constrains(std::unordered_map<unsigned short, std::vector<Eigen::Vector3i>>)
                    {
                        {
                            0,
                            {
                                {0, 4, 7},
                                {1, 5, 6},
                                {5, 6, 7},
                                {4, 5, 7},
                                {1, 0, 4},
                                {1, 5, 4},
                                {1, 0, 6},
                                {0, 6, 7},
                            }
                        },
                        {
                            1,
                            {
                                {1, 2, 6},
                                {0, 3, 7},
                                {2, 6, 7},
                                {2, 3, 7},
                                {1, 0, 2},
                                {0, 2, 3},
                                {1, 0, 6},
                                {0, 6, 7},
                            }
                        }
                    },
                }
            },
            // case 6.1.2
            {
                // connected vertices(unordered_set<std::unordered_set<Vertex>>)
                {{Vertex::v1, Vertex::v6}, {Vertex::v2, Vertex::v4}},
                // FeatureMC33Table
                {
                    // feature interpolation rules(std::vector<std::tuple<unsigned short, unsigned short>>)
                    {},
                    // feature triangles(std::unordered_map<unsigned short, std::vector<unsigned short>>)
                    {},
                    // mc33 triangles(std::vector<unsigned short>)
                    {},
                    // common triangles(std::vector<unsigned short>)
                    {0x59A, 0xA91, 0x198, 0x182, 0x287, 0x27A, 0xA75},
                    // feature point constrains(std::unordered_map<unsigned short, std::vector<Eigen::Vector3i>>)
                    {},
                }
            },
            // case 6.2
            {
                // connected vertices(unordered_set<std::unordered_set<Vertex>>)
                {{Vertex::v2, Vertex::v4}, {Vertex::v2, Vertex::v5}},
                // FeatureMC33Table
                {
                    // feature interpolation rules(std::vector<std::tuple<unsigned short, unsigned short>>)
                    {{5, 0x21987}},
                    // feature triangles(std::unordered_map<unsigned short, std::vector<unsigned short>>)
                    {
                        {0, {0x21, 0x19, 0x98, 0x87, 0x72}},
                    },
                    // mc33 triangles(std::vector<unsigned short>)
                    {},
                    // common triangles(std::vector<unsigned short>)
                    {},
                    // feature point constrains(std::unordered_map<unsigned short, std::vector<Eigen::Vector3i>>)
                    {
                        {
                            0,
                            CUBE_TRIANGLES
                        }
                    },
                }
            },
        }
    },

    // case 7.1
    {
        // distance signs
        0b01011011,
        // cube
        {
            // case 7.1
            {
                // connected vertices(unordered_set<std::unordered_set<Vertex>>)
                {{Vertex::v3, Vertex::v6}, {Vertex::v4, Vertex::v6}, {Vertex::v1, Vertex::v6}, {Vertex::v0, Vertex::v6}},
                // FeatureMC33Table
                {
                    // feature interpolation rules(std::vector<std::tuple<unsigned short, unsigned short>>)
                    {{3, 0x12A}, {3, 0x67B}, {3, 0x459}},
                    // feature triangles(std::unordered_map<unsigned short, std::vector<unsigned short>>)
                    {
                        {0, {0x21, 0xA2, 0x1A}},
                        {1, {0x67, 0x7B, 0xB6}},
                        {2, {0x45, 0x59, 0x94}}
                    },
                    // mc33 triangles(std::vector<unsigned short>)
                    {},
                    // common triangles(std::vector<unsigned short>)
                    {},
                    // feature point constrains(std::unordered_map<unsigned short, std::vector<Eigen::Vector3i>>)
                    {
                        {
                            0,
                            {
                                {2, 3, 6},
                                {1, 2, 6},
                                {0, 1, 2},
                                {0, 2, 3},
                                {0, 3, 6},
                                {0, 1, 6},
                            }
                        },
                        {
                            1,
                            {
                                {3, 6, 7},
                                {4, 6, 7},
                                {0, 3, 7},
                                {0, 4, 7},
                                {0, 3, 6},
                                {0, 4, 6},
                            }
                        },
                        {
                            2,
                            {
                                {4, 5, 6},
                                {1, 5, 6},
                                {0, 1, 4},
                                {1, 4, 5},
                                {0, 6, 4},
                                {0, 6, 1},
                            }
                        }
                    },
                }
            },
            // case 7.2
            {
                // connected vertices(unordered_set<std::unordered_set<Vertex>>)
                {{Vertex::v2, Vertex::v7}, {Vertex::v4, Vertex::v6}, {Vertex::v1, Vertex::v6}, {Vertex::v0, Vertex::v6}},
                // FeatureMC33Table
                {
                    // feature interpolation rules(std::vector<std::tuple<unsigned short, unsigned short>>)
                    {{6, 0xB76A12}, {3, 0x459}},
                    // feature triangles(std::unordered_map<unsigned short, std::vector<unsigned short>>)
                    {
                        {0, {0xB7, 0x76, 0x6A, 0xA1, 0x12, 0x2B}},
                        {1, {0x45, 0x59, 0x94}}
                    },
                    // mc33 triangles(std::vector<unsigned short>)
                    {},
                    // common triangles(std::vector<unsigned short>)
                    {},
                    // feature point constrains(std::unordered_map<unsigned short, std::vector<Eigen::Vector3i>>)
                    {
                        {
                            0,
                            {
                                {4, 7, 6},
                                {1, 2, 6},
                                {1, 0, 4},
                                {0, 3, 7},
                                {0, 7, 4},
                                {0, 1, 2},
                                {0, 2, 3},
                                {2, 3, 6},
                                {3, 6, 7},
                                {1, 4, 6},
                            }
                        },
                        {
                            1,
                            {
                                {4, 5, 6},
                                {1, 5, 6},
                                {1, 4, 5},
                                {1, 4, 6},
                            }
                        },
                    },
                }
            },
            // case 7.3
            {
                // connected vertices(unordered_set<std::unordered_set<Vertex>>)
                {{Vertex::v2, Vertex::v7}, {Vertex::v2, Vertex::v5}, {Vertex::v4, Vertex::v6}, {Vertex::v0, Vertex::v6}},
                // FeatureMC33Table
                {
                    // feature interpolation rules(std::vector<std::tuple<unsigned short, unsigned short>>)
                    {},
                    // feature triangles(std::unordered_map<unsigned short, std::vector<unsigned short>>)
                    {},
                    // mc33 triangles(std::vector<unsigned short>)
                    {0x7B, 0xB2, 0x21, 0x19, 0x94, 0x45, 0x5A, 0xA6, 0x67},
                    // common triangles(std::vector<unsigned short>)
                    {},
                    // feature point constrains(std::unordered_map<unsigned short, std::vector<Eigen::Vector3i>>)
                    {},
                }
            },
            // case 7.4.1
            {
                // connected vertices(unordered_set<std::unordered_set<Vertex>>)
                {{Vertex::v2, Vertex::v7}, {Vertex::v2, Vertex::v5}, {Vertex::v5, Vertex::v7}},
                // FeatureMC33Table
                {
                    // feature interpolation rules(std::vector<std::tuple<unsigned short, unsigned short>>)
                    {{3, 0x56A}, {6, 0x2B7491}},
                    // feature triangles(std::unordered_map<unsigned short, std::vector<unsigned short>>)
                    {
                        {0, {0x65, 0xA6, 0x5A}},
                        {1, {0xB2, 0x7B, 0x47, 0x94, 0x19, 0x21}}
                    },
                    // mc33 triangles(std::vector<unsigned short>)
                    {},
                    // common triangles(std::vector<unsigned short>)
                    {},
                    // feature point constrains(std::unordered_map<unsigned short, std::vector<Eigen::Vector3i>>)
                    {
                        {
                            0,
                            {
                                {2, 7, 6},
                                {5, 6, 7},
                                {2, 6, 5},
                                {2, 5, 7},
                            }
                        },
                        {
                            1,
                            {
                                {2, 3, 7},
                                {4, 5, 7},
                                {1, 2, 5},
                                {0, 3, 7},
                                {0, 4, 7},
                                {0, 1, 3},
                                {1, 2, 3},
                                {0, 1, 4},
                                {1, 4, 5},
                                {2, 5, 7},
                            }
                        },
                    },
                }
            },
            // case 7.4.2
            {
                // connected vertices(unordered_set<std::unordered_set<Vertex>>)
                {{Vertex::v2, Vertex::v7}, {Vertex::v2, Vertex::v5}, {Vertex::v5, Vertex::v7}, {Vertex::v0, Vertex::v6}},
                // FeatureMC33Table
                {
                    // feature interpolation rules(std::vector<std::tuple<unsigned short, unsigned short>>)
                    {},
                    // feature triangles(std::unordered_map<unsigned short, std::vector<unsigned short>>)
                    {},
                    // mc33 triangles(std::vector<unsigned short>)
                    {},
                    // common triangles(std::vector<unsigned short>)
                    {0x67B, 0x657, 0x547, 0x459, 0x519, 0x5A1, 0x1A2, 0x2A6, 0x26B},
                    // feature point constrains(std::unordered_map<unsigned short, std::vector<Eigen::Vector3i>>)
                    {},
                }
            },
        }
    },

    // case 8
    {
        // distance signs
        0b11001100,
        // cube
        {
            {
                // connected vertices(unordered_set<std::unordered_set<Vertex>>)
                {},
                // FeatureMC33Table
                {
                    // feature interpolation rules(std::vector<std::tuple<unsigned short, unsigned short>>)
                    {{4, 0x1573}},
                    // feature triangles(std::unordered_map<unsigned short, std::vector<unsigned short>>)
                    {
                        {0, {0x51, 0x75, 0x37, 0x13}}
                    },
                    // mc33 triangles(std::vector<unsigned short>)
                    {},
                    // common triangles(std::vector<unsigned short>)
                    {},
                    // feature point constrains(std::unordered_map<unsigned short, std::vector<Eigen::Vector3i>>)
                    {
                        {
                            0,
                            CUBE_TRIANGLES
                        }
                    },
                }
            }
        }
    },

    // case 9
    {
        // distance signs
        0b11100100,
        // cube
        {
            {
                // connected vertices(unordered_set<std::unordered_set<Vertex>>)
                {},
                // FeatureMC33Table
                {
                    // feature interpolation rules(std::vector<std::tuple<unsigned short, unsigned short>>)
                    {{6, 0xB74912}},
                    // feature triangles(std::unordered_map<unsigned short, std::vector<unsigned short>>)
                    {
                        {0, {0xB7, 0x74, 0x49, 0x91, 0x12, 0x2B}}
                    },
                    // mc33 triangles(std::vector<unsigned short>)
                    {},
                    // common triangles(std::vector<unsigned short>)
                    {},
                    // feature point constrains(std::unordered_map<unsigned short, std::vector<Eigen::Vector3i>>)
                    {
                        {
                            0,
                            CUBE_TRIANGLES
                        }
                    },
                }
            }
        }
    },

    // case 10
    {
        // distance signs
        0b01101001,
        // cube
        {
            // 10.1.1
            {
                // connected vertices(unordered_set<std::unordered_set<Vertex>>)
                {{Vertex::v3, Vertex::v5}, {Vertex::v0, Vertex::v5}, {Vertex::v0, Vertex::v6}, {Vertex::v3, Vertex::v5}},
                // FeatureMC33Table
                {
                    // feature interpolation rules(std::vector<std::tuple<unsigned short, unsigned short>>)
                    {{4, 0xB648}, {4, 0x2A90}},
                    // feature triangles(std::unordered_map<unsigned short, std::vector<unsigned short>>)
                    {
                        {0, {0xB6, 0x64, 0x48, 0x8B}},
                        {1, {0xA2, 0x9A, 0x09, 0x20}}
                    },
                    // mc33 triangles(std::vector<unsigned short>)
                    {},
                    // common triangles(std::vector<unsigned short>)
                    {},
                    // feature point constrains(std::unordered_map<unsigned short, std::vector<Eigen::Vector3i>>)
                    {
                        {
                            0,
                            {
                                {3, 6, 7},
                                {0, 5, 4},
                                FACE_TRIANGLES.at(Face::f3)[0],
                                FACE_TRIANGLES.at(Face::f3)[1],
                                FACE_TRIANGLES.at(Face::f5)[0],
                                FACE_TRIANGLES.at(Face::f5)[1],
                            }
                        },
                        {
                            1,
                            {
                                {2, 3, 6},
                                {0, 1, 5},
                                FACE_TRIANGLES.at(Face::f1)[0],
                                FACE_TRIANGLES.at(Face::f1)[1],
                                FACE_TRIANGLES.at(Face::f4)[0],
                                FACE_TRIANGLES.at(Face::f4)[1],
                            }
                        }
                    },
                }
            },
            // 10.1.2
            {
                // connected vertices(unordered_set<std::unordered_set<Vertex>>)
                {{Vertex::v3, Vertex::v5}, {Vertex::v0, Vertex::v2}, {Vertex::v2, Vertex::v4}, {Vertex::v1, Vertex::v7}},
                // FeatureMC33Table
                {
                    // feature interpolation rules(std::vector<std::tuple<unsigned short, unsigned short>>)
                    {},
                    // feature triangles(std::unordered_map<unsigned short, std::vector<unsigned short>>)
                    {},
                    // mc33 triangles(std::vector<unsigned short>)
                    {},
                    // common triangles(std::vector<unsigned short>)
                    {0xB6A, 0x2BA, 0x20B, 0xB08, 0x048, 0x409, 0x64A, 0x49A},
                    // feature point constrains(std::unordered_map<unsigned short, std::vector<Eigen::Vector3i>>)
                    {},
                }
            },
            // 10.2
            {
                // connected vertices(unordered_set<std::unordered_set<Vertex>>)
                {{Vertex::v0, Vertex::v5}, {Vertex::v1, Vertex::v4}, {Vertex::v2, Vertex::v4}, {Vertex::v1, Vertex::v7}},
                // FeatureMC33Table
                {
                    // feature interpolation rules(std::vector<std::tuple<unsigned short, unsigned short>>)
                    {},
                    // feature triangles(std::unordered_map<unsigned short, std::vector<unsigned short>>)
                    {},
                    // mc33 triangles(std::vector<unsigned short>)
                    {0xA2, 0x9A, 0x49, 0x64, 0xB6, 0x8B, 0x08},
                    // common triangles(std::vector<unsigned short>)
                    {},
                    // feature point constrains(std::unordered_map<unsigned short, std::vector<Eigen::Vector3i>>)
                    {},
                }
            }
        }
    },

    // case 11
    {
        // distance signs
        0b11101000,
        // cube
        {
            {
                // connected vertices(unordered_set<std::unordered_set<Vertex>>)
                {},
                // FeatureMC33Table
                {
                    // feature interpolation rules(std::vector<std::tuple<unsigned short, unsigned short>>)
                    {{6, 0x23749A}},
                    // feature triangles(std::unordered_map<unsigned short, std::vector<unsigned short>>)
                    {
                        {0, {0x23, 0x37, 0x74, 0x49, 0x9A, 0xA2}}
                    },
                    // mc33 triangles(std::vector<unsigned short>)
                    {},
                    // common triangles(std::vector<unsigned short>)
                    {},
                    // feature point constrains(std::unordered_map<unsigned short, std::vector<Eigen::Vector3i>>)
                    {
                        {
                            0,
                            CUBE_TRIANGLES
                        }
                    },
                }
            }
        }
    },

    // case 12
    {
        // distance signs
        0b01011100,
        // cube
        {
            // 12.1.1
            {
                // connected vertices(unordered_set<std::unordered_set<Vertex>>)
                {{Vertex::v3, Vertex::v4}, {Vertex::v4, Vertex::v6}},
                // FeatureMC33Table
                {
                    // feature interpolation rules(std::vector<std::tuple<unsigned short, unsigned short>>)
                    {{3, 0x67B}, {5, 0x13845}},
                    // feature triangles(std::unordered_map<unsigned short, std::vector<unsigned short>>)
                    {
                        {0, {0x67, 0x7B, 0xB6}},
                        {1, {0x13, 0x38, 0x84, 0x45, 0x51}},
                    },
                    // mc33 triangles(std::vector<unsigned short>)
                    {},
                    // common triangles(std::vector<unsigned short>)
                    {},
                    // feature point constrains(std::unordered_map<unsigned short, std::vector<Eigen::Vector3i>>)
                    {
                        {
                            0,
                            {
                                {3, 6, 7},
                                {4, 6, 7},
                                {3, 4, 7},
                                {3, 4, 6}
                            }
                        },
                        {
                            1,
                            {
                                {2, 3, 6},
                                {4, 5, 6},
                                {0, 3, 4},
                                FACE_TRIANGLES.at(Face::f1)[0],
                                FACE_TRIANGLES.at(Face::f1)[1],
                                FACE_TRIANGLES.at(Face::f2)[0],
                                FACE_TRIANGLES.at(Face::f2)[1],
                                FACE_TRIANGLES.at(Face::f4)[0],
                                FACE_TRIANGLES.at(Face::f4)[1],
                                {3, 4, 6}
                            }
                        }
                    },
                }
            },
            // 12.1.2
            {
                // connected vertices(unordered_set<std::unordered_set<Vertex>>)
                {{Vertex::v3, Vertex::v4}, {Vertex::v4, Vertex::v6}, {Vertex::v1, Vertex::v7}},
                // FeatureMC33Table
                {
                    // feature interpolation rules(std::vector<std::tuple<unsigned short, unsigned short>>)
                    {},
                    // feature triangles(std::unordered_map<unsigned short, std::vector<unsigned short>>)
                    {},
                    // mc33 triangles(std::vector<unsigned short>)
                    {},
                    // common triangles(std::vector<unsigned short>)
                    {0x38B, 0xB87, 0x478, 0x467, 0x456, 0x516, 0xB61, 0xB13},
                    // feature point constrains(std::unordered_map<unsigned short, std::vector<Eigen::Vector3i>>)
                    {},
                }
            },
            // 12.2
            {
                // connected vertices(unordered_set<std::unordered_set<Vertex>>)
                {{Vertex::v3, Vertex::v4}, {Vertex::v5, Vertex::v7}, {Vertex::v1, Vertex::v7}},
                // FeatureMC33Table
                {
                    // feature interpolation rules(std::vector<std::tuple<unsigned short, unsigned short>>)
                    {},
                    // feature triangles(std::unordered_map<unsigned short, std::vector<unsigned short>>)
                    {},
                    // mc33 triangles(std::vector<unsigned short>)
                    {0x51, 0x13, 0x38, 0x84, 0x47, 0x7B, 0xB6, 0x65},
                    // common triangles(std::vector<unsigned short>)
                    {},
                    // feature point constrains(std::unordered_map<unsigned short, std::vector<Eigen::Vector3i>>)
                    {},
                }
            },
            // 12.3
            {
                // connected vertices(unordered_set<std::unordered_set<Vertex>>)
                {{Vertex::v4, Vertex::v6}, {Vertex::v2, Vertex::v4}, {Vertex::v0, Vertex::v7}},
                // FeatureMC33Table
                {
                    // feature interpolation rules(std::vector<std::tuple<unsigned short, unsigned short>>)
                    {},
                    // feature triangles(std::unordered_map<unsigned short, std::vector<unsigned short>>)
                    {},
                    // mc33 triangles(std::vector<unsigned short>)
                    {0x84, 0x45, 0x51, 0x13, 0x3B, 0xB6, 0x67, 0x78},
                    // common triangles(std::vector<unsigned short>)
                    {},
                    // feature point constrains(std::unordered_map<unsigned short, std::vector<Eigen::Vector3i>>)
                    {},
                }
            },
        }
    },

    // case 13
    {
        // distance signs
        0b10100101,
        // cube
        {
            // case 13.1
            {
                // connected vertices(unordered_set<std::unordered_set<Vertex>>)
                {{Vertex::v2, Vertex::v7}, {Vertex::v0, Vertex::v5}, {Vertex::v0, Vertex::v7}, {Vertex::v2, Vertex::v5}, {Vertex::v5, Vertex::v7}, {Vertex::v0, Vertex::v2}},
                // FeatureMC33Table
                {
                    // feature interpolation rules(std::vector<std::tuple<unsigned short, unsigned short>>)
                    {{3, 0x23B}, {3, 0x56A}, {3, 0x478}, {3, 0x019}},
                    // feature triangles(std::unordered_map<unsigned short, std::vector<unsigned short>>)
                    {
                        {0, {0x32, 0xB3, 0x2B}},
                        {1, {0x56, 0x6A, 0xA5}},
                        {2, {0x48, 0x87, 0x74}},
                        {3, {0x10, 0x91, 0x09}},
                    },
                    // mc33 triangles(std::vector<unsigned short>)
                    {},
                    // common triangles(std::vector<unsigned short>)
                    {},
                    // feature point constrains(std::unordered_map<unsigned short, std::vector<Eigen::Vector3i>>)
                    {
                        {
                            0,
                            {
                                {0, 2, 3},
                                {0, 3, 7},
                                {2, 3, 7},
                                {static_cast<int>(Vertex::vc), 2, 7},
                                {static_cast<int>(Vertex::vc), 0, 7},
                                {static_cast<int>(Vertex::vc), 0, 2},
                            }
                        },
                        {
                            1,
                            {
                                {2, 6, 7},
                                {2, 6, 5},
                                {5, 6, 7},
                                {static_cast<int>(Vertex::vc), 2, 7},
                                {static_cast<int>(Vertex::vc), 7, 5},
                                {static_cast<int>(Vertex::vc), 5, 2},
                            }
                        },
                        {
                            2,
                            {
                                {0, 4, 7},
                                {0, 4, 5},
                                {5, 4, 7},
                                {static_cast<int>(Vertex::vc), 0, 5},
                                {static_cast<int>(Vertex::vc), 5, 7},
                                {static_cast<int>(Vertex::vc), 0, 7},
                            }
                        },
                        {
                            3,
                            {
                                {0, 1, 2},
                                {5, 1, 2},
                                {0, 1, 5},
                                {static_cast<int>(Vertex::vc), 0, 2},
                                {static_cast<int>(Vertex::vc), 2, 5},
                                {static_cast<int>(Vertex::vc), 0, 5},
                            }
                        }
                    },
                }
            },
            // case 13.2
            {
                // connected vertices(unordered_set<std::unordered_set<Vertex>>)
                {{Vertex::v3, Vertex::v6}, {Vertex::v0, Vertex::v5}, {Vertex::v0, Vertex::v7}, {Vertex::v2, Vertex::v5}, {Vertex::v5, Vertex::v7}, {Vertex::v0, Vertex::v2}},
                // FeatureMC33Table
                {
                    // feature interpolation rules(std::vector<std::tuple<unsigned short, unsigned short>>)
                    {{6, 0x23B65A}, {3, 0x478}, {3, 0x019}},
                    // feature triangles(std::unordered_map<unsigned short, std::vector<unsigned short>>)
                    {
                        {0, {0x32, 0xB3, 0x6B, 0x56, 0xA5, 0x2A}},
                        {1, {0x48, 0x87, 0x74}},
                        {2, {0x10, 0x91, 0x09}},
                    },
                    // mc33 triangles(std::vector<unsigned short>)
                    {},
                    // common triangles(std::vector<unsigned short>)
                    {},
                    // feature point constrains(std::unordered_map<unsigned short, std::vector<Eigen::Vector3i>>)
                    {
                        {
                            1,
                            {
                                {0, 4, 7},
                                {0, 4, 5},
                                {4, 5, 7},
                                {0, 7, 5},
                            }
                        },
                        {
                            2,
                            {
                                {0, 1, 2},
                                {5, 1, 2},
                                {0, 1, 5},
                                {0, 2, 5},
                            }
                        },
                        {
                            0,
                            {
                                {5, 6, 7},
                                {2, 6, 5},
                                {0, 3, 7},
                                {0, 2, 3},
                                FACE_TRIANGLES.at(Face::f2)[0],
                                FACE_TRIANGLES.at(Face::f2)[1],
                                {0, 7, 5},
                                {0, 2, 5},
                            }
                        }
                    },
                }
            },
            // case 13.3
            {
                // connected vertices(unordered_set<std::unordered_set<Vertex>>)
                {{Vertex::v3, Vertex::v6}, {Vertex::v0, Vertex::v5}, {Vertex::v3, Vertex::v4}, {Vertex::v2, Vertex::v5}, {Vertex::v5, Vertex::v7}, {Vertex::v0, Vertex::v2}},
                // FeatureMC33Table
                {
                    // feature interpolation rules(std::vector<std::tuple<unsigned short, unsigned short>>)
                    {{3, 0x019}},
                    // feature triangles(std::unordered_map<unsigned short, std::vector<unsigned short>>)
                    {
                        {0, {0x10, 0x91, 0x09}},
                    },
                    // mc33 triangles(std::vector<unsigned short>)
                    {0xA5, 0x2A, 0x32, 0x83, 0x48, 0x74, 0xB7, 0x6B, 0x56},
                    // common triangles(std::vector<unsigned short>)
                    {},
                    // feature point constrains(std::unordered_map<unsigned short, std::vector<Eigen::Vector3i>>)
                    {
                        {
                            0,
                            {
                                {2, 1, 5},
                                {2, 1, 0},
                                {0, 1, 5},
                                {2, 0, 5},
                            }
                        }
                    },
                }
            },
            // case 13.4
            {
                // connected vertices(unordered_set<std::unordered_set<Vertex>>)
                {{Vertex::v0, Vertex::v6}, {Vertex::v0, Vertex::v5}, {Vertex::v3, Vertex::v4}, {Vertex::v2, Vertex::v5}, {Vertex::v5, Vertex::v7}, {Vertex::v1, Vertex::v3}},
                // FeatureMC33Table
                {
                    // feature interpolation rules(std::vector<std::tuple<unsigned short, unsigned short>>)
                    {},
                    // feature triangles(std::unordered_map<unsigned short, std::vector<unsigned short>>)
                    {},
                    // mc33 triangles(std::vector<unsigned short>)
                    {0x30, 0x83, 0x48, 0x74, 0xB7, 0x6B, 0x56, 0xA5, 0x2A, 0x12, 0x91, 0x09},
                    // common triangles(std::vector<unsigned short>)
                    {},
                    // feature point constrains(std::unordered_map<unsigned short, std::vector<Eigen::Vector3i>>)
                    {},
                }
            },
            // case 13.5.1
            {
                // connected vertices(unordered_set<std::unordered_set<Vertex>>)
                {{Vertex::v3, Vertex::v6}, {Vertex::v0, Vertex::v5}, {Vertex::v0, Vertex::v7}, {Vertex::v1, Vertex::v6}, {Vertex::v5, Vertex::v7}, {Vertex::v1, Vertex::v3}},
                // FeatureMC33Table
                {
                    // feature interpolation rules(std::vector<std::tuple<unsigned short, unsigned short>>)
                    {{3, 0x12A}, {3, 0x478}, {6, 0x3B6590}},
                    // feature triangles(std::unordered_map<unsigned short, std::vector<unsigned short>>)
                    {
                        {0, {0x12, 0x2A, 0xA1}},
                        {1, {0x48, 0x87, 0x74}},
                        {2, {0x3B, 0xB6, 0x65, 0x59, 0x90, 0x03}}
                    },
                    // mc33 triangles(std::vector<unsigned short>)
                    {},
                    // common triangles(std::vector<unsigned short>)
                    {},
                    // feature point constrains(std::unordered_map<unsigned short, std::vector<Eigen::Vector3i>>)
                    {
                        {
                            0,
                            {
                                {1, 2, 3},
                                {1, 2, 6},
                                {6, 2, 3},
                                {1, 3, 6},
                            }
                        },
                        {
                            1,
                            {
                                {0, 4, 7},
                                {0, 4, 5},
                                {5, 4, 7},
                                {0, 7, 5},
                            }
                        },
                        {
                            2,
                            {
                                {
                                    {0, 7, 3},
                                    {0, 3, 1},
                                    {0, 1, 5},
                                    {1, 5, 6},
                                    {5, 6, 7},
                                    {3, 6, 7},
                                    {0, 7, 5},
                                    {1, 3, 6},
                                }
                            }
                        }
                    },
                }
            },
            // case 13.5.2
            {
                // connected vertices(unordered_set<std::unordered_set<Vertex>>)
                {{Vertex::v3, Vertex::v5}, {Vertex::v1, Vertex::v4}, {Vertex::v0, Vertex::v7}, {Vertex::v1, Vertex::v6}, {Vertex::v5, Vertex::v7}, {Vertex::v1, Vertex::v3}},
                // FeatureMC33Table
                {
                    // feature interpolation rules(std::vector<std::tuple<unsigned short, unsigned short>>)
                    {{3, 0x12A}},
                    // feature triangles(std::unordered_map<unsigned short, std::vector<unsigned short>>)
                    {
                        {0, {0x12, 0x2A, 0xA1}},
                    },
                    // mc33 triangles(std::vector<unsigned short>)
                    {},
                    // common triangles(std::vector<unsigned short>)
                    {0x480, 0x830, 0x387, 0x37B, 0x76B, 0x674, 0x645, 0x495, 0x094},
                    // feature point constrains(std::unordered_map<unsigned short, std::vector<Eigen::Vector3i>>)
                    {
                        {
                            0,
                            {
                                {1, 2, 3},
                                {1, 2, 6},
                                {6, 2, 3},
                                {1, 3, 6},
                            }
                        },
                    },
                }
            }
        }
    },

    // case 14
    {
        // distance signs
        0b11010100,
        // cube
        {
            {
                // connected vertices(unordered_set<std::unordered_set<Vertex>>)
                {},
                // FeatureMC33Table
                {
                    // feature interpolation rules(std::vector<std::tuple<unsigned short, unsigned short>>)
                    {{6, 0x2B8451}},
                    // feature triangles(std::unordered_map<unsigned short, std::vector<unsigned short>>)
                    {
                        {0, {0x2B, 0xB8, 0x84, 0x45, 0x51, 0x12}}
                    },
                    // mc33 triangles(std::vector<unsigned short>)
                    {},
                    // common triangles(std::vector<unsigned short>)
                    {},
                    // feature point constrains(std::unordered_map<unsigned short, std::vector<Eigen::Vector3i>>)
                    {
                        {
                            0,
                            CUBE_TRIANGLES
                        }
                    },
                }
            }
        }
    },
};

std::ostream& operator<<(std::ostream& stream, const std::unordered_set<std::unordered_set<Vertex>, boost::hash<std::unordered_set<Vertex>>>& object)
{
    for(const auto& connected_vertice: object)
    {
        for(const auto& connected_vertex: connected_vertice)
        {
            stream << static_cast<int>(connected_vertex) << ",";
        }
        stream << std::endl;
    }
    return stream;
}

TEST(GlobalTest, init_tables)
{
    init_tables();
    print_mc33_table();
}

TEST(GlobalTest, vertex_interpolation_connected_cube_diagnal_connected)
{
    std::vector<double> signed_distance1{-100, 1, 1, 1, 1, 1, -100, 1};
    Vertex v1 = Vertex::v0;
    Vertex v2 = Vertex::v6;
    auto connected1 = vertex_interpolation_connected(signed_distance1, v1, v2);
    ASSERT_EQ(connected1, true);
}

TEST(GlobalTest, vertex_interpolation_connected_cube_diagnal_disconnected)
{
    std::vector<double> signed_distance1{-1, 100, 100, 100, 100, 100, -1, 100};
    Vertex v1 = Vertex::v0;
    Vertex v2 = Vertex::v6;
    auto connected1 = vertex_interpolation_connected(signed_distance1, v1, v2);
    ASSERT_EQ(connected1, false);
}

TEST(GlobalTest, vertex_interpolation_connected_face_diagnal_connected)
{
    std::vector<double> signed_distance1{-100, 1, -100, 101, 100, 100, -1, 100};
    Vertex v1 = Vertex::v0;
    Vertex v2 = Vertex::v2;
    auto connected1 = vertex_interpolation_connected(signed_distance1, v1, v2);
    ASSERT_EQ(connected1, true);
}

TEST(GlobalTest, vertex_interpolation_connected_face_diagnal_disconnected)
{
    std::vector<double> signed_distance1{-1, 100, -1, 100, 100, 100, -1, 100};
    Vertex v1 = Vertex::v0;
    Vertex v2 = Vertex::v2;
    auto connected1 = vertex_interpolation_connected(signed_distance1, v1, v2);
    ASSERT_EQ(connected1, false);
}

TEST(GlobalTest, vertex_connected_same_edge_connected)
{
    std::bitset<8> signed_distance{0b11111100};
    Vertex v1 = Vertex::v0;
    Vertex v2 = Vertex::v1;
    auto connected = vertex_connected(signed_distance, v1, v2);
    ASSERT_EQ(connected, true);
}

TEST(GlobalTest, vertex_connected_same_edge_disconnected)
{
    std::bitset<8> signed_distance{0b11111110};
    Vertex v1 = Vertex::v0;
    Vertex v2 = Vertex::v1;
    auto connected = vertex_connected(signed_distance, v1, v2);
    ASSERT_EQ(connected, false);
}

TEST(GlobalTest, vertex_connected_edge_diagnal_connected)
{
    std::bitset<8> signed_distance{0b11111000};
    Vertex v1 = Vertex::v0;
    Vertex v2 = Vertex::v2;
    auto connected = vertex_connected(signed_distance, v1, v2);
    ASSERT_EQ(connected, true);
}

TEST(GlobalTest, vertex_connected_edge_diagnal_disconnected)
{
    {
        std::bitset<8> signed_distance{0b11111010};
        Vertex v1 = Vertex::v0;
        Vertex v2 = Vertex::v2;
        auto connected = vertex_connected(signed_distance, v1, v2);
        ASSERT_EQ(connected, false);
    }

    {
        std::bitset<8> signed_distance{0b11001011};
        Vertex v1 = Vertex::v1;
        Vertex v2 = Vertex::v6;
        auto connected = vertex_connected(signed_distance, v1, v2);
        ASSERT_EQ(connected, false);
    }
}

TEST(GlobalTest, vertex_connected_cube_diagnal_connected)
{
    std::bitset<8> signed_distance{0b10111000};
    Vertex v1 = Vertex::v0;
    Vertex v2 = Vertex::v2;
    auto connected = vertex_connected(signed_distance, v1, v2);
    ASSERT_EQ(connected, true);
}

TEST(GlobalTest, vertex_connected_cube_diagnal_disconnected)
{
    {
        std::bitset<8> signed_distance{0b01111000};
        Vertex v1 = Vertex::v1;
        Vertex v2 = Vertex::v7;
        auto connected = vertex_connected(signed_distance, v1, v2);
        ASSERT_EQ(connected, false);
    }

    {
        std::bitset<8> signed_distance{0b10111110};
        Vertex v1 = Vertex::v0;
        Vertex v2 = Vertex::v6;
        auto connected = vertex_connected(signed_distance, v1, v2);
        ASSERT_EQ(connected, false);
    }
}

TEST(GlobalTest, has_value_bigger_than_zero_in_interval)
{
    double start = 0.0f, end = 1.0f;
    // test a bigger than 0, there is no value than 0
    {
        double a = 1.0f, b = 1.0f, c = -3.0f;
        auto result = has_value_bigger_than_zero_in_interval(a, b, c, start, end);
        ASSERT_EQ(result, false);
    }
    // test a bigger than 0, there are values bigger than 0
    // start is bigger than 0, but end not
    {
        double a = 1.0f, b = -3.0f, c = 1.0f;
        auto result = has_value_bigger_than_zero_in_interval(a, b, c, start, end);
        ASSERT_EQ(result, true);
    }
    // end is bigger than 0, but start not
    {
        double a = 1.0f, b = 1.0f, c = -1.0f;
        auto result = has_value_bigger_than_zero_in_interval(a, b, c, start, end);
        ASSERT_EQ(result, true);
    }
    // start and end are all bigger than 0
    {
        double a = 1.0f, b = 1.0f, c = 1.0f;
        auto result = has_value_bigger_than_zero_in_interval(a, b, c, start, end);
        ASSERT_EQ(result, true);
    }

    // test a smaller than 0
    // test no root
    {
        double a = -1.0f, b = -1.0f, c = -1.0f;
        auto result = has_value_bigger_than_zero_in_interval(a, b, c, start, end);
        ASSERT_EQ(result, false);
    }
    // test has root, but no value bigger than 0
    {
        double a = -1.0f, b = -3.0f, c = -1.0f;
        auto result = has_value_bigger_than_zero_in_interval(a, b, c, start, end);
        ASSERT_EQ(result, false);
    }
    // test has root, and has value bigger than 0
    {
        double a = -1.0f, b = 1.0f, c = -1.0f / 8.0f;
        auto result = has_value_bigger_than_zero_in_interval(a, b, c, start, end);
        ASSERT_EQ(result, true);
    }
}

TEST(GlobalTest, feature_point_triangle_same_length)
{
    for(const auto&[distances_sign, sub_cases]: MC33_TABLES)
    {
        for(const auto& [connected_vertices, feature_mc33_table]: sub_cases)
        {
            EXPECT_EQ(feature_mc33_table.feature_interpolation_rules.size(), feature_mc33_table.fp_connected_edges.size()) <<  "distances_sign: " << distances_sign << ", connected_vertices: " << connected_vertices << " feature mc33 table has problem in feature length equality";
            EXPECT_EQ(feature_mc33_table.feature_interpolation_rules.size(), feature_mc33_table.constrains.size()) <<  "distances_sign: " << distances_sign << ", connected_vertices: " << connected_vertices << " feature mc33 table has problem in feature length equality";
        }
    }
}