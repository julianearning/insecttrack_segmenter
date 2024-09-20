#ifndef ASSIGN_HPP
#define ASSIGN_HPP

#include <vector>
#include "../utils/pointcloud.hpp"

class Assign {
private:
public:
    Assign(PointCloud* downsampled, PointCloud* full, std::vector<int>* labels_downsampled, std::vector<int>* labels_full, size_t ndims);
    Assign(PointCloud* downsampled, PointCloud* full, std::vector<int>* labels_downsampled, std::vector<int>* labels_full, std::vector<size_t>* noise_points, size_t ndims);
};

#endif