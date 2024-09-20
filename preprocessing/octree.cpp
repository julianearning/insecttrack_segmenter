#include "octree.hpp"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <execution>
#include <limits>

#define DIMENSIONX 1280
#define DIMENSIONY 720

OcTree::OcTree(PointCloud * pc, PointCloud * output, int target_size, ssize_t min_points, size_t ndims) {
    this->pc = pc;
    for(size_t i = 0; i<pc->size(); i++) {
        nodes.push_back(OcTreeNode(pc->at(i), i));
    }
    this->output = output;
    this->target_size = target_size;
    this->min_points = min_points;
    this->ndims = ndims;
}


void OcTree::partitiont() {
    int last_pivot=0;
    pivot_indexes_t.push_back(last_pivot);
    for(int i = 1; i<(int)(pc->size()); i++) {
        if((pc->at(i).t - pc->at(last_pivot).t) > target_size) {
            pivot_indexes_t.push_back(i);
            last_pivot=i;
        }
    }
    pivot_indexes_t.push_back(pc->size());
} 



bool compareBySpace(const OcTreeNode& p1, const OcTreeNode& p2, const size_t dim) {
    p1.p.data.at(dim) < p2.p.data.at(dim);
}


bool compareByT(const Point& p1, const Point& p2) { 
    return p1.t < p2.t; 
}

bool compareByIdx(const OcTreeNode& p1, const OcTreeNode& p2) {
    return p1.idx < p2.idx;
}

//Point OcTree::centroid(PointCloud * c) {
Point OcTree::centroid(std::vector<OcTreeNode>* ocn) {
    Point p=Point(ndims);
    size_t n = ocn->size();
    for(int i = 0; i<(int)n; i++) {
        p = p + ocn->at(i).p;
    }
    p = p / n;
    return p;
}


//Point OcTree::closest_to_centroid(PointCloud * c) {
Point OcTree::closest_to_centroid(std::vector<OcTreeNode>* ocn) {
    Point p=Point(ndims);
    size_t n = ocn->size();
    for(int i = 0; i<(int)n; i++) {
        p = p + ocn->at(i).p;
    }
    p = p / n;
    Point result=Point(ndims);
    double min_distance=std::numeric_limits<double>::max();
    double curr_distance=0.0;
    for(int i = 0; i<(int)ocn->size(); i++) {
        curr_distance = ocn->at(i).p.no_root_euclidian_distance(p);
        if(curr_distance<min_distance) {
            min_distance=curr_distance;
            result=ocn->at(i).p;
        }
    }
    return result;
}


void OcTree::partition_space(std::vector<OcTreeNode>* buffer1, size_t dim) {
    Point p;
    int last_pivot;
    std::vector<OcTreeNode> buffer2;
    std::sort( std::execution::par, buffer1->begin(), buffer1->end(), [dim](const OcTreeNode& p1, const OcTreeNode& p2) {
                  return p1.p.data.at(dim) < p2.p.data.at(dim);
              } );
    last_pivot=0;
    for(int j = 0; j<((int)buffer1->size()); j++) {
        if((buffer1->at(j).p.data.at(dim)-buffer1->at(last_pivot).p.data.at(dim)) > target_size) {
            last_pivot=j;
            if(dim == (ndims-2)) {
                if((buffer2.size()>min_points)) {
                    p = closest_to_centroid(&buffer2);
                    output->push_back(p);
                } else {
                    std::sort(buffer2.begin(), buffer2.end(), compareByIdx);
                    for(size_t i = 0; i<buffer2.size(); i++) {
                        noise_points.push_back(buffer2.at(i).idx);
                    }
                }
            } else {
                partition_space(&buffer2, dim+1);
            }
            buffer2.clear();
        }
        buffer2.push_back(buffer1->at(j));
        if(j==((int)buffer1->size()-1)) {  // Rest
            if(dim == (ndims-2)) {
                if((buffer2.size()>min_points)) {
                    p = closest_to_centroid(&buffer2);
                    output->push_back(p);
                } else {
                    std::sort(buffer2.begin(), buffer2.end(), compareByIdx);
                    for(size_t i = 0; i<buffer2.size(); i++) {
                        noise_points.push_back(buffer2.at(i).idx);
                    }
                }
            } else {
                partition_space(&buffer2, dim+1);
            }
            buffer2.clear();
        }
    }
}

void OcTree::partition_first_space_dimension() {
    int first_pt;
    int last_pt;
    int last_pivot;
    size_t dim=0;
    //PointCloud slice_copy;
    //PointCloud buffer1;
    std::vector<OcTreeNode> slice_copy;
    std::vector<OcTreeNode> buffer1;
    std::vector<int> pivot_indexes_x;
    Point p=Point(ndims);
    for(int t = 0; t<((int)pivot_indexes_t.size()-1); t++) { // for time index 
        first_pt = pivot_indexes_t.at(t);
        last_pt = pivot_indexes_t.at(t+1)-1;
        slice_copy.clear();
        pivot_indexes_x.clear();
        for(int i = first_pt; i<=last_pt; i++) {
            slice_copy.push_back(nodes.at(i));
        }
        std::stable_sort( std::execution::par, slice_copy.begin(), slice_copy.end(), [dim](const OcTreeNode& p1, const OcTreeNode& p2) {
                  return p1.p.data.at(dim) < p2.p.data.at(dim);
              });
        last_pivot=0;
        buffer1.clear();
        for(int i = 0; i<(int)slice_copy.size(); i++) {
            if(((slice_copy.at(i).p.data.at(0) - slice_copy.at(last_pivot).p.data.at(0)) > target_size) || i==((int)slice_copy.size()-1)) {
                last_pivot=i;

                partition_space(&buffer1, 1);

                buffer1.clear();
            } 
            buffer1.push_back(slice_copy.at(i));
        } 
        if(buffer1.size()>0) {
            partition_space(&buffer1, 1);
            buffer1.clear();
        }
    }
}  

void OcTree::preprocess(std::vector<size_t>* denoised_points) {
    partitiont();
    //partitionx();
    partition_first_space_dimension();
    std::sort(std::execution::par, output->begin(), output->end(), compareByT);
    std::sort(noise_points.begin(), noise_points.end());
    for(size_t i = 0; i<noise_points.size(); i++) {
        denoised_points->push_back(noise_points.at(i));
    }
}
