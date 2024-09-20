#ifndef MST_HPP
#define MST_HPP

#include "../utils/pointcloud.hpp"
#include <vector>



class UndirectedNode {
public:
    std::vector<UndirectedNode *> neighbours;
    std::vector<float> distances;
    size_t index;
    UndirectedNode(size_t index) {this->index=index;};
};


class UndirectedMST {
private:
    ssize_t num_clusters;
    std::vector<UndirectedNode> mst;
    PointCloud* pc;
    double max_mean_edge_len=-1;
    size_t min_cluster_size;
    size_t min_depth;
    size_t n;
    float max_edge_dist;
    float max_branch_depth;
    std::vector<int> * labels;
    void remove_long_edges(double max_edge_len);
    void remove_local_long_edges(size_t depth, double c);
    void edge_dist_at_depth(UndirectedNode * node, UndirectedNode * taboo_node, size_t depth, double* mean, double* sd);
    float get_edge_len(std::pair<size_t, size_t> edge);
    size_t color_nodes();
    double expand(UndirectedNode * node, int cluster_id); // returns mean edge dist
    //void merge_labels(size_t label_a, size_t label_b);
    void remove_branches();
    bool probe_depth(size_t node_index, size_t node_neighbour, size_t max_depth);
    void reconnect(size_t num_pca_points);
public:
    UndirectedMST(PointCloud* pc, std::vector<int> * labels, float max_edge_dist, size_t min_cluster_size, const size_t ndims, double max_mean_edge_len);  // used for clustering
    UndirectedMST(PointCloud* pc, std::vector<int> * labels, size_t min_depth, size_t min_cluster_size, const size_t ndims);   // used for postprocessing
    size_t colorMST();  // returns the number of labels generated, input is true if branches should be removed. When false long edges are removed.
    bool draw_mst(const char * fname, bool is4d);
};



class DirectedNode {
public:
    std::vector<DirectedNode *> predecessor;
    DirectedNode * successor;
    std::vector<double> dist_predecessor;
    double dist_successor;
    int index;
    DirectedNode(int index) {this->successor=NULL; this->index=index;  };
};


class DirectedMST {
private:
    int octree_target_size;
    PointCloud removed_edges;
    double max_edge_dist;
    std::vector<DirectedNode> mst;
    std::vector<int> * labels;
    size_t min_cluster_size;
    PointCloud * pc;
    void build_MST();
    ssize_t BFS(size_t index, PointCloud * n1, PointCloud * n2, size_t max_depth);   // undirected breadth-first search for both sides of edge mst.at(index) to mst.at(index)->successor, returns Pointcloud for both sides (before edge=n1, after edge=n2)
    //void color_nodes();
    int expand(DirectedNode * node, int cluster_id);
    //double calc_max_edge_dist_estimate(std::vector<std::vector<float>> * dists);
    void calc_pca(std::vector<std::vector<size_t>> * dists);
    double inconsistent_edge_length(size_t index_a , size_t index_b, std::vector<float> * dists_a, std::vector<float> * dists_b, std::vector<size_t> * indices_a, std::vector<size_t> * indices_b);
    void remove_long_edges();
    double choose_global_max_dist_threshold();
    double choose_local_max_dist_threshold(size_t index);
    void delete_edge(size_t index);
public:
    DirectedMST(PointCloud * pc, std::vector<int> * labels, double max_edge_dist, int octree_target_size, size_t min_cluster_size);
    bool draw_mst(const char * fname, bool is4d);
    size_t colorMST();
}; 


#endif