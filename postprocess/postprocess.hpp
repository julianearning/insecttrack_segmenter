#ifndef POSTPROCESS_HPP
#define POSTPROCESS_HPP

#include "../utils/pointcloud.hpp"
#include <vector>

class PostProcessor {
private:
    PointCloud* pc;
    std::vector<int>* labels;
    float threshold_sd;
    float min_depth;
    size_t max_label;
    float point_proportion;
    double tolerance;
    double min_cosine_similarity;
    size_t min_cluster_size;
    size_t ndims;
    void get_cluster(size_t c, PointCloud* input, std::vector<int>* input_labels, PointCloud* cluster);
    float compute_spacial_variance(PointCloud* cluster);
    void get_pca_of_cluster(std::pair<PointCloud,PointCloud>* segments, Point0 * a_begin, Point0* b_begin, Point0* a_end, Point0* b_end);
    //void get_velocity_of_cluster(std::pair<PointCloud,PointCloud>* segments, Point0 *a_end, Point0 *b_end); 
public:
    PostProcessor(PointCloud * pc, std::vector<int> * labels, float threshold_sd, float min_depth, size_t max_label, float point_proportion, double tolerance, double min_cosine_similarity, size_t min_cluster_size, size_t ndims);
    size_t postprocess();
};

#endif 