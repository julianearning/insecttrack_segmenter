#ifndef MATHSTUFF_HPP
#define MATHSTUFF_HPP

#include "pointcloud.hpp"

class MathStuff {
public:
    double orthogonal_lsq(PointCloud & pc, Point* a, Point* b, size_t ndims);
    double orthogonal_lsq_distance(Point* a, Point* b, Point * c, size_t ndims);
};

#endif