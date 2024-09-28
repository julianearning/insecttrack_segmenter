#ifndef OCTREE_HPP
#define OCTREE_HPP

#include "../utils/pointcloud.hpp"


class OcTreeNode {
public:
    Point0 p;
    size_t idx;
    OcTreeNode() = default;
    OcTreeNode(Point0 p, size_t idx) { this->p=p; this->idx=idx; }
    OcTreeNode(const OcTreeNode& ocn) { this->p=ocn.p; this->idx=ocn.idx; }
    OcTreeNode(size_t ndims) { p = Point0(ndims); idx=0; }
};


class OcTree {
private:
    std::vector<OcTreeNode> nodes;
    PointCloud * pc;
    PointCloud * output;
    std::vector<size_t> noise_points;
    int dimensiont;  // 
    int dimensionx;// 
    int dimensiony; // 
    int target_size;
    size_t ndims;
    ssize_t min_points;
    void partitiont();
    //void partitionx();
    void partition();
    //void partitiony(std::vector<OcTreeNode>* buffer);
    void partition_space(std::vector<OcTreeNode>* buffer, size_t dim);
    void partition_first_space_dimension();
    std::vector<int> pivot_indexes_t;   // every element is the index of a nex box in t direction
    std::vector<int> * indexes_left;
    //int binary_search_t(int start_idx, int end_idx, double search_t);
    //Point0 centroid(PointCloud * pc);
    //Point0 closest_to_centroid(PointCloud * pc);
    Point0 centroid(std::vector<OcTreeNode>* ocn);
    Point0 closest_to_centroid(std::vector<OcTreeNode>* ocn);
public:
    OcTree(PointCloud * pc, PointCloud * output, int target_size, ssize_t min_points, size_t ndims);
    void preprocess(std::vector<size_t> * denoised_points);
};


#endif