#ifndef MATHSTUFF_HPP
#define MATHSTUFF_HPP

#include "pointcloud.hpp"

class MathStuff {
public:
    double orthogonal_lsq(PointCloud & pc, Point0* a, Point0* b, size_t ndims);
    double orthogonal_lsq_distance(Point0* a, Point0* b, Point0 * c, size_t ndims);
};

#endif