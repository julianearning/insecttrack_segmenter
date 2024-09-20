#ifndef FILEIO_HPP
#define FILEIO_HPP

#include <string>
#include <vector>
#include "pointcloud.hpp"

class FileParser {
    private:
        std::vector<double> original_time_information;
        size_t find_pos_start(std::ifstream & is, Point & p, bool is4d);
    public:
        int parse_file(std::string filename, PointCloud & pc, double scale_factor_t, bool is4d);    // read in file and output pc
        int write_csv(std::string filename, PointCloud & pc, std::vector<int> & labels, bool append, std::vector<size_t>* noise_points);
        int write_csv(std::string filename, PointCloud & pc, std::vector<int> & labels, bool append);
        int write_loess(std::string filename, std::vector<PointCloud> * data, std::vector<int> * labels, bool append, double scale_factor_t);
        int write_csv_no_label(std::string filename, PointCloud & pc);
        int merge_csv(std::string filename, PointCloud & pc, std::vector<int> & labels, double scale_factor_t, bool is4d);   // for window mode, TODO: implementation
        bool debug_gnuplot(PointCloud &cloud, PointCloud &cloud_downsampled, const char *fname, bool is4d); 
        //int clusters_to_gnuplot(const PointCloud &cloud, const std::vector<int> &labels);
};


#endif