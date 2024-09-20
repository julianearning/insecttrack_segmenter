#ifndef OCTREE_HPP
#define OCTREE_HPP

#include "pointcloud.hpp"

class BoundingBox {
public:
    Point upper_left;
    int dimension;
    BoundingBox& operator=(const BoundingBox& other);
    BoundingBox(){};
    BoundingBox(Point upper_left, int dimension) {this->upper_left=upper_left; this->dimension=dimension;};
    BoundingBox(int upper_left_x, int upper_left_y, double upper_left_t, int dimension) { this->upper_left.x=upper_left_x; this->upper_left.y=upper_left.y; this->upper_left.t=upper_left.t; this->dimension=dimension;};
};


class OcTree {
private:
    BoundingBox m_region;
    PointCloud m_objects;
    PointCloud pending;
    OcTree * m_childNode0;
    OcTree * m_childNode1;
    OcTree * m_childNode2;
    OcTree * m_childNode3;
    OcTree * m_childNode4;
    OcTree * m_childNode5;
    OcTree * m_childNode6;
    OcTree * m_childNode7;
    char m_activeNodes=0;
    int max_depth;
    int curr_depth;
    OcTree * _parent;
public:
    OcTree(PointCloud & pc, int max_depth, int curr_depth, BoundingBox & bb);
    //OcTree(Point upper_left, int dimension, int max_depth, int curr_depth);
};


#endif