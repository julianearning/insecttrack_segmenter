#ifndef MST_HPP
#define MST_HPP

#include "../utils/pointcloud.hpp"
#include <vector>

class DirectedNode {
public:
    std::vector<DirectedNode *> predecessor;
    DirectedNode * successor;
    std::vector<double> dist_predecessor;
    double dist_successor;
    int index;
    DirectedNode(int index) {this->successor=NULL; this->index=index;};
};


class DirectedMST {
private:
    int octree_target_size;
    PointCloud removed_edges;
    double max_edge_dist;
    double pca_angle_value;
    size_t max_branch_depth;
    std::vector<DirectedNode> mst;
    std::vector<int> * labels;
    std::vector<std::vector<Point>> pca_equations;
    std::vector<std::vector<Point>> pca_equations_n1;
    std::vector<std::vector<Point>> pca_equations_n2;
    std::vector<double> kneeness;
    PointCloud * pc;
    bool pca_angle_criterion;
    bool edge_length_criterion;
    bool twoD;
    void build_MST();
    //void build_MST_twoD();
    bool DFS(DirectedNode * n, size_t curr_depth);   // directed depth-first search until depth, recursive implementation 
    ssize_t BFS(size_t index, PointCloud * n1, PointCloud * n2, size_t max_depth);   // undirected breadth-first search for both sides of edge mst.at(index) to mst.at(index)->successor, returns Pointcloud for both sides (before edge=n1, after edge=n2)
    void color_nodes();
    int expand(DirectedNode * node, int cluster_id);
    double calc_max_edge_dist_estimate(std::vector<std::vector<float>> * dists);
    void calc_pca(std::vector<std::vector<size_t>> * dists);
    double inconsistent_edge_length(size_t index_a , size_t index_b, std::vector<float> * dists_a, std::vector<float> * dists_b, std::vector<size_t> * indices_a, std::vector<size_t> * indices_b);
    void remove_long_edges();
    void remove_time_gaps();
    void remove_knees();           // knees: Knick in der Bahn
    void remove_branches();
    double choose_global_max_dist_threshold();
    double choose_local_max_dist_threshold(size_t index);
    void delete_edge(size_t index);
public:
    DirectedMST(PointCloud * pc, std::vector<int> * labels, double max_edge_dist, size_t max_branch_depth, bool edge_length_criterion, bool pca_angle_criterion, double pca_angle_value, int octree_target_size);
    bool draw_mst(const char * fname);
    bool draw_pca(const char * fname);
    bool draw_pca_debug(const char * fname, bool twoD);
};


class UndirectedNode {
public:
    std::vector<UndirectedNode *> neighbours;
    std::vector<double> dist_heighbours;
    int index;
    UndirectedNode(int index) {this->index=index;};
};


struct Edge {
    size_t i;
    size_t j;
    float dist;
};


class UndirectedMST {
private:
    int octree_target_size;
    PointCloud removed_edges;
    double max_edge_dist;
    double pca_angle_value;
    size_t max_branch_depth;
    std::vector<UndirectedNode> mst;
    std::vector<int> * labels;
    PointCloud * pc;
    bool pca_angle_criterion;
    bool edge_length_criterion;
public:
    UndirectedMST(PointCloud * pc, std::vector<int> * labels, double max_edge_dist, size_t max_branch_depth, bool edge_length_criterion, bool pca_angle_criterion, double pca_angle_value, int octree_target_size);
    bool draw_mst(const char * fname);
    void build_MST();
};



#endif