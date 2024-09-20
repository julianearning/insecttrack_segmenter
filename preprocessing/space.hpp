#ifndef SPACE_HPP
#define SPACE_HPP

#include "../utils/pointcloud.hpp"

class Space {
private:
    PointCloud * pc;
    PointCloud * output;
    double eps;
    bool pseudorandom_order=false;
public:
    Space(PointCloud * pc, PointCloud * output, double eps);
    void preprocess();
};

#endif