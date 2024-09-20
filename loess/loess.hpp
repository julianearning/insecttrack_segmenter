#ifndef LOESS_HPP
#define LOESS_HPP

#include "../utils/pointcloud.hpp"
#include <vector>

class LOESS {
    size_t ndims;
public:
    LOESS(PointCloud * pc, std::vector<int> * labels, std::vector<PointCloud> * loess, std::vector<int> * label_order, double span, size_t ndims);
    void do_loess(PointCloud * pc, PointCloud * out, double span);
};


#endif