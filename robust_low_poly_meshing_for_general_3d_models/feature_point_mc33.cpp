#include <vector>
#include <iterator>
#include <numeric>
#include <fstream>
#include <unordered_set>
#include <bitset>
#include <chrono>

#define eigen_assert(X) do { if(!(X)) throw std::runtime_error(#X); } while(false);
#include "Eigen/Dense"
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include "igl/readOFF.h"
#include "igl/signed_distance.h"
#include "igl/AABB.h"
#include "igl/per_face_normals.h"

#include "feature_point_mc33.h"
#include "feature_point_mc33_tables.h"
#include "unit_test.h"
#include "charles_mc33_type.h"
#include "multivariable_extream.h"
#include "feature_mc33_cache.h"

//#include <boost/stacktrace.hpp>

// interpolation cache
FeatureMC33Cache<Eigen::Vector3d> FEATURE_MC33_CACHE;
int cache_hitted_count = 0;

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
    int z, int y, int x,
    Edge edge,
    const std::vector<Eigen::Vector3d>& coors,
    const std::vector<double>& signed_distance
)
{
    try{
        auto cache_value = FEATURE_MC33_CACHE.get({z, y, x, edge});
        std::cout << "cache hitted" << std::endl;
        cache_hitted_count++;
        return cache_value;
    }catch(const std::out_of_range &exception)
    {
        // do nothing
    }

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
    FEATURE_MC33_CACHE.set({z, y, x, edge}, vertex);
    return vertex;
}

FeatureMC33Table  get_feature_mc33_table(const std::vector<double>& signed_distances)
{
    std::bitset<8> distances_signs{0b11111111};
    for(int i = 0; i < signed_distances.size(); i++)
    {
        auto signed_distance = signed_distances[i];
        if(signed_distance < 0)
        {
            distances_signs.set(i, false);
        }
    }
    auto sub_cases = MC33_TABLES.at(distances_signs);
    // get all connected vertices that are not ambiguous
    std::unordered_set<std::unordered_set<Vertex>, boost::hash<std::unordered_set<Vertex>>> connected_verticess;
    for(unsigned short i = 0; i < 8; i++)
    {
        for(unsigned short j = i + 1; j < 8; j++)
        {
            Vertex vertex1 = static_cast<Vertex>(i), vertex2 = static_cast<Vertex>(j);
            if (distances_signs[static_cast<unsigned short>(vertex1)] == distances_signs[static_cast<unsigned short>(vertex2)])
            {
                if (!vertex_connected(distances_signs, vertex1, vertex2))
                {
                    if (vertex_interpolation_connected(signed_distances, vertex1, vertex2))
                    {
                        connected_verticess.emplace(std::unordered_set<Vertex>{vertex1, vertex2});
                    }
                }
            }
        }
    }

    // get egess by connected vertices
    auto feature_mc33_table = sub_cases.at(connected_verticess);
    return feature_mc33_table;
}

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

class Chaos
{
public:
    // vertex of source mesh
    Eigen::MatrixXd V;
    // face of source mesh
    Eigen::MatrixXi F;
    // normal of face of source mesh
    Eigen::MatrixXd FN;
    // AABB tree of source mesh
    igl::AABB<Eigen::MatrixXd, 3> tree;

    void init_source_mesh(const std::string& off_path)
    {
        igl::readOFF(off_path, this->V, this->F);
        this->tree.init(this->V, this->F);
        igl::per_face_normals(this->V, this->F, FN);
    }

    Eigen::Vector3d get_interpolated_vertex_normal(const Eigen::Vector3d& interpolated_vertex)
    {
        Eigen::Vector3d normal;

        Eigen::VectorXd sqrD;
        Eigen::VectorXi I;
        Eigen::MatrixXd C;
        Eigen::MatrixXd P(1, 3);
        P(0, 0) = interpolated_vertex[0];
        P(0, 1) = interpolated_vertex[1];
        P(0, 2) = interpolated_vertex[2];

        //for(int i = 0; i < P.rows(); i++)
        //{
        //    for(int j = 0; j < P.cols(); j++)
        //    {
        //        P(i, j) = interpolated_vertex(i * 2 + j);
        //    }
        //}
        // igl::point_mesh_squared_distance(P,V,F,sqrD,I,C);
        tree.squared_distance(this->V, this->F, P, sqrD, I, C);

        if(sqrD[0] == 0)
        {
            // use normal of triangle of mesh directly
            normal[0] = this->FN.row(I[0])[0];
            normal[1] = this->FN.row(I[0])[1];
            normal[2] = this->FN.row(I[0])[2];
        }
        else
        {
            // normal = interpolated_vertex - C.row(I[0]);
            normal[0] = interpolated_vertex[0] - C.row(0)[0];
            normal[1] = interpolated_vertex[1] - C.row(0)[1];
            normal[2] = interpolated_vertex[2] - C.row(0)[2];
        }

        return normal;
    }

    Eigen::Vector3d get_feature_vertex(
        std::vector<Eigen::Vector3d>& interpolated_vertices, std::vector<Eigen::Vector3d>& normals,
        const std::vector<Eigen::Vector3d>& vertices, const std::vector<Eigen::Vector3i>& triangles
    )
    {
        Eigen::Vector3d feature_vertex;

        // get expression
        auto x = SymEngine::symbol("x");
        auto y = SymEngine::symbol("y");
        auto z = SymEngine::symbol("z");
        SymEngine::Expression x_(x);
        SymEngine::Expression y_(y);
        SymEngine::Expression z_(z);
        SymEngine::Expression ex;

        for(auto iter = interpolated_vertices.begin(); iter != interpolated_vertices.end(); iter++)
        {
            auto index = std::distance(interpolated_vertices.begin(), iter);
            auto interpolated_point = *iter;
            auto interpolated_normal = normals[index];
            ex = ex + SymEngine::pow(interpolated_normal.x() * (x_ - interpolated_point.x()) + interpolated_normal.y() * (y_ - interpolated_point.y()) + interpolated_normal.z() * (z_ - interpolated_point.z()), 2);
        }

        auto [min_value, max_value, min_location, max_location] = min_max_3_variable_quadratic_polynomial(ex, {x, y, z}, vertices, triangles);

        feature_vertex[0] = min_location.at(x);
        feature_vertex[1] = min_location.at(y);
        feature_vertex[2] = min_location.at(z);

        return feature_vertex;
    }

    Eigen::Vector3d get_feature_vertex(
        int z, int y, int x,
        const std::tuple<unsigned short, unsigned short>& feature_interpolation_rule,
        const std::vector<Eigen::Vector3d>& coors,
        const std::vector<double>& signed_distances,
        const std::vector<Eigen::Vector3i>& feature_constrain
    )
    {
        auto [length, edges] = feature_interpolation_rule;
        std::vector<Eigen::Vector3d> interpolated_vertices;
        std::vector<Eigen::Vector3d> normals;
        for(int i = 0; i < length; i++)
        {
            edges = edges >> (i == 0 ? 0 : 4);
            Edge edge = static_cast<Edge>(edges & 0xF);
            auto vertex = get_vertex_by_edge(z, y, x, edge, coors, signed_distances);
            auto normal = get_interpolated_vertex_normal(vertex);
            interpolated_vertices.emplace_back(std::move(vertex));
            normals.emplace_back(normal);
        }
        // get feature point from normals and interpolated vertices
        Eigen::Vector3d center_point = (coors[0] + coors[6]) / 2;
        std::vector<Eigen::Vector3d> vertices;
        vertices.insert(vertices.end(), coors.begin(), coors.end());
        vertices.emplace_back(std::move(center_point));
        auto feature_vertex = get_feature_vertex(interpolated_vertices, normals, vertices, feature_constrain);

        // may be additional test is needed, the result is abnormal, some point seems out of the cube
        // NOTE: test result: passed, no point outside the cube, maybe my concern is redunant
        // if(feature_vertex.x() < coors[0].x() || feature_vertex.x() > coors[6].x() || feature_vertex.y() < coors[0].y() || feature_vertex.y() > coors[6].y() || feature_vertex.z() < coors[0].z() || feature_vertex.z() > coors[6].z())
        // {
        //     throw std::exception("feature vertex out of cube! please take attention!");
        // }

        return feature_vertex;
    }

    std::tuple<std::vector<Eigen::Vector3d>, std::vector<Eigen::Vector3i>> get_feature_mc33_triangles(
        int z, int y, int x,
        const FeatureMC33Table& feature_mc33_table,
        const std::vector<Eigen::Vector3d>& coors,
        const std::vector<double>& signed_distances
    )
    {
        std::vector<Eigen::Vector3d> vertices;
        std::vector<Eigen::Vector3i> triangles;

        std::vector<Eigen::Vector3d> feature_vertices;
        // for(const auto& feature_interpolation_rule: feature_mc33_table.feature_interpolation_rules)
        for(auto iter = feature_mc33_table.feature_interpolation_rules.begin(); iter != feature_mc33_table.feature_interpolation_rules.end(); iter++)
        {
            auto feature_interpolation_rule = *iter;
            auto index = std::distance(feature_mc33_table.feature_interpolation_rules.begin(), iter);
            auto feature_constrain = feature_mc33_table.constrains.at(index);
            auto feature_vertice = get_feature_vertex(z, y, x, feature_interpolation_rule, coors, signed_distances, feature_constrain);
            feature_vertices.emplace_back(feature_vertice);
        }
        auto lambda_func = [&](const Edge& edge) -> Eigen::Vector3d {
            auto vertex = get_vertex_by_edge(z, y, x, edge, coors, signed_distances);
            return vertex;
        };
        std::function<Eigen::Vector3d(const Edge&)> func = lambda_func;
        // compute feature triangles
        for(const auto& [feature_index, f_triangles]: feature_mc33_table.fp_connected_edges)
        {
            for(const auto& f_triangle: f_triangles)
            {
                auto temp_vertices = handle_edges(f_triangle, FeatureMC33Table::feature_triangle_length, func);
                temp_vertices.emplace_back(feature_vertices[feature_index]);
                if(temp_vertices.size() != 3)
                {
                    throw std::exception("triangle that contains feature vertex doesn't has three vertice");
                }
                int current_vertex_index = vertices.size();
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
            int current_vertex_index = vertices.size();
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
            int current_vertex_index = vertices.size();
            triangles.emplace_back(Eigen::Vector3i{current_vertex_index, current_vertex_index + 1, current_vertex_index + 2});
            vertices.insert(vertices.end(), temp_vertices.begin(), temp_vertices.end());
        }
        return {vertices, triangles};
    }

    void run()
    {
        int nx = 5;
        int ny = 5;
        int nz = 5;

        Eigen::Vector3d m = V.colwise().minCoeff();
        Eigen::Vector3d M = V.colwise().maxCoeff();

        Eigen::MatrixXd P;
        P.resize(nx * ny * nz, 3);
        for(unsigned int k = 0; k < nz; k++)
        {
            double z_axis = m.z() + double(k) * (M.z() - m.z()) / double(nz - 1);
            for(unsigned int j = 0; j < ny; j++)
            {
                double y_axis = m.y() + double(j) * (M.y() - m.y()) / double(ny - 1);
                for(unsigned int i = 0; i < nx; i++)
                {
                    double x_axis = m.x() + double(i) * (M.x() - m.x()) / double(nx - 1);
                    unsigned int count = k * ny * nx + j * nx + i;
                    P(count, 0) = x_axis;
                    P(count, 1) = y_axis;
                    P(count, 2) = z_axis;
                }
            }
        }

        Eigen::VectorXi I;
        Eigen::MatrixXd N,C;
        Eigen::VectorXd S;
        igl::SignedDistanceType type = igl::SIGNED_DISTANCE_TYPE_PSEUDONORMAL;
        igl::signed_distance(P, V, F, type, S, I, C, N);

        // std::cout << "------------------S---------------------" << std::endl;
        // std::cout << S << std::endl;

        std::vector<Eigen::Vector3d> vertices;
        std::vector<Eigen::Vector3i> triangles;
        int iso_value = 0;

        for(int k = 0; k < nz - 1; k++)
        {
            for(int j = 0; j < ny - 1; j++)
            {
                for(int i = 0; i < nx - 1; i++)
                {
                    //  get cube
                    std::vector<Eigen::Vector3d> coors{
                        {static_cast<double>(i), static_cast<double>(j), static_cast<double>(k)},
                        {static_cast<double>(i), static_cast<double>(j + 1), static_cast<double>(k)},
                        {static_cast<double>(i), static_cast<double>(j + 1), static_cast<double>(k + 1)},
                        {static_cast<double>(i), static_cast<double>(j), static_cast<double>(k + 1)},
                        {static_cast<double>(i + 1), static_cast<double>(j), static_cast<double>(k)},
                        {static_cast<double>(i + 1), static_cast<double>(j + 1), static_cast<double>(k)},
                        {static_cast<double>(i + 1), static_cast<double>(j + 1), static_cast<double>(k + 1)},
                        {static_cast<double>(i + 1), static_cast<double>(j), static_cast<double>(k + 1)},
                    };
                    unsigned int count = k * ny * nx + j * nx + i;
                    std::vector<double> signed_distance{
                        S[k * ny * nx + j * nx + i],
                        S[k * ny * nx + (j + 1) * nx + i],
                        S[(k + 1) * ny * nx + (j + 1) * nx + i],
                        S[(k + 1) * ny * nx + j * nx + i],
                        S[k * ny * nx + j * nx + i + 1],
                        S[k * ny * nx + (j + 1) * nx + i + 1],
                        S[(k + 1) * ny * nx + (j + 1) * nx + i + 1],
                        S[(k + 1) * ny * nx + j * nx + i + 1],
                    };

                    FeatureMC33Table feature_mc33_table;
                    try
                    {
                        // edgess = MC_TABLES.at(distance_symbol);
                        feature_mc33_table = get_feature_mc33_table(signed_distance);
                    }
                    catch (const std::out_of_range& e)
                    {
                        std::cerr << "Caught std::out_of_range exception: " << e.what() << std::endl;
                        throw e;
                    }
                    auto [temp_vertices, temp_triangles] = this->get_feature_mc33_triangles(k, j, i, feature_mc33_table, coors, signed_distance);
                    for(int l = 0; l < temp_triangles.size(); l++)
                    {
                        auto current_vertices_size = vertices.size();
                        temp_triangles[l][0] = temp_triangles[l][0] + current_vertices_size;
                        temp_triangles[l][1] = temp_triangles[l][1] + current_vertices_size;
                        temp_triangles[l][2] = temp_triangles[l][2] + current_vertices_size;
                    }
                    vertices.insert(vertices.end(), temp_vertices.begin(), temp_vertices.end());
                    triangles.insert(triangles.end(), temp_triangles.begin(), temp_triangles.end());
                }
            }
        }
        write_obj("./bunny.obj", vertices, triangles);
    }
};


TEST(GlobalTest, bunny)
{
    init_tables();
    Chaos chaos;
    chaos.init_source_mesh("D:\\Library\\libigl\\build\\_deps\\libigl_tutorial_tata-src\\bunny.off");
    auto start = std::chrono::high_resolution_clock::now();
    chaos.run();
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "total cost seconds is " << duration.count() << "s" << std::endl;
    std::cout << "cache hitted count" << cache_hitted_count << std::endl;
    std::cout << "cache count" << FEATURE_MC33_CACHE.size() << std::endl;
}