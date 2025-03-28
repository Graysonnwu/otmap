#ifndef MESH_H
#define MESH_H

#include <cmath>
#include <vector>
#include <string>
#include <numeric>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>

struct HashPair {
    template <class T1, class T2>
    std::size_t operator () (const std::pair<T1, T2>& p) const {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);

        // Using XOR (^) to combine hashes
        return h1 ^ h2;
    }
};

class Mesh {
    private:

        void generate_structured_mesh(int nx, int ny, double width, double height, std::vector<std::vector<int>> &triangles, std::vector<std::vector<double>> &points);
        void build_vertex_to_triangles();

    public:
        Mesh(double width, double height, int res_x, int res_y);
        ~Mesh();

        std::unordered_map<int, std::vector<int>> vertex_to_triangles;

        std::vector<std::vector<double>> source_points;
        std::vector<std::vector<double>> target_points;

        std::vector<std::vector<int>> triangles;

        double width;
        double height;

        int res_x;
        int res_y;

        std::pair<std::vector<std::pair<int, int>>, std::vector<int>> find_adjacent_elements(int vertex_index);

        void find_vertex_connectivity(int vertex_index, std::vector<int> &neighborList, std::vector<int> &neighborMap, const std::vector<int>& mask_indices);

        std::vector<std::vector<double>> calculate_refractive_normals(double focal_len, double refractive_index);
        std::vector<std::vector<double>> calculate_refractive_normals_uniform(std::vector<std::vector<double>> target_pts, double focal_len, double refractive_index); 

        void save_solid_obj_target(double thickness, const std::string& filename);
        void save_solid_obj_source(double thickness, const std::string& filename);

        std::vector<std::vector<double>> circular_transform(std::vector<std::vector<double>> &input_points);

        void get_vertex_neighbor_ids(int vertex_id, int &left_vertex, int &right_vertex, int &top_vertex, int &bottom_vertex);

        std::vector<double> calculate_vertex_normal(std::vector<std::vector<double>> &points, int vertex_index);

        bool is_border(int vertex_id);
};

#endif
