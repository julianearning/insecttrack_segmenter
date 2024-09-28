#ifndef DBSCAN_HPP
#define DBSCAN_HPP

#include "../utils/pointcloud.hpp"
#include <vector>

class DBSCAN {
    private:
        PointCloud * pc_ptr; 
        int limited_calc_DBSCAN(std::vector<int> * labels_out, const double eps, const int minPts, const int buffer_size, const double max_lsq_distance);
        int vanilla_calc_DBSCAN(std::vector<int> * labels_out, const double eps, const int minPts);
    public:
        DBSCAN(PointCloud * pc_ptr);
        int calc_DBSCAN(std::vector<int> * labels_out, const double eps, const int minPts, const bool vanilla, const bool limited, const int buffer_size, const double max_lsq_distance);
        //double orthogonal_lsq(PointCloud & pc, Point0* a, Point0* b);
};


#endif
