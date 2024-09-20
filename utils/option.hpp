#ifndef OPTION_HPP
#define OPTION_HPP
#include <string>

class Opt {
    private:
        double eps;
        int minPts;
        double space_eps;
        double max_lsq_distance;
        double min_cosine_score;
        int buffer_size;
        double scale;
        bool smooth;
        bool remove_duplicates;
        double smooth_radius;
        int mst_max_branch_depth;
        double mst_max_edge_dist;
        bool mst_pca_angle;
        double mst_pca_angle_value;
        bool mst_twoD;
        size_t mst_min_cluster_size;
        double mst_max_mean_edge_len;
        std::string preprocessing_method;
        std::string segmentation_method;
        std::string out_filename;
        int octree_target_size;
        int octree_denoise;
        std::string octree_selection;
        bool generate_preprocessing_debug;
        std::string generate_preprocessing_debug_filename;
        bool generate_mst_debug_gnuplot;
        std::string mst_debug_gnuplot_filename;
        bool window_mode;
        bool silent;
        bool no_postprocessing;
        float postprocessing_threshold_sd;
        //ssize_t postprocessing_min_depth;
        float postprocessing_min_depth;
        float postprocessing_point_proportion;
        double postprocessing_tolerance;
        double postprocessing_min_cosine_similarity;
        bool is4d;
    public:
        Opt();
        int parse_args(int argc, char ** argv);
        double get_scale() { return this->scale; };
        double get_max_lsq_distance() { return this->max_lsq_distance; };
        bool get_smooth() { return this->smooth; };
        std::string get_preprocessing_method() {return this->preprocessing_method;};
        int get_octree_target_size() { return this->octree_target_size; }; 
        std::string get_segmentation_method() { return this->segmentation_method;};
        std::string get_octree_selection() { return this->octree_selection;};
        double get_dbscan_eps() { return this->eps; };
        int get_dbscan_minPts() {return this->minPts; }; 
        int get_buffer_size() { return this->buffer_size; };
        int get_space_eps() { return this->space_eps; };
        int get_smooth_radius() { return this->smooth_radius; };
        double get_min_cosine_score() { return this->min_cosine_score; }; 
        double get_mst_max_edge_dist() { return this->mst_max_edge_dist; };
        int get_mst_max_branch_depth() { return this->mst_max_branch_depth; };
        bool get_generate_preprocessing_debug() { return this->generate_preprocessing_debug; };
        std::string get_generate_preprocessing_debug_filename() { return this->generate_preprocessing_debug_filename; };
        std::string get_out_filename() { return this->out_filename; };
        bool get_generate_mst_debug_gnuplot() { return this->generate_mst_debug_gnuplot; };
        std::string get_mst_debug_gnuplot_filename() { return this->mst_debug_gnuplot_filename; };
        bool get_remove_duplicates() { return this->remove_duplicates; };
        bool get_window_mode() { return this->window_mode; };
        bool get_mst_pca_angle() { return this->mst_pca_angle; };
        double get_mst_pca_angle_value() { return this->mst_pca_angle_value; };
        bool get_mst_twoD() { return this->mst_twoD; }; 
        int get_octree_denoise() { return this->octree_denoise; };
        bool get_silent() { return this->silent; };
        float get_postprocessing_threshold_sd() { return this->postprocessing_threshold_sd; };
        float get_postprocessing_point_proportion() { return this->postprocessing_point_proportion; }; 
        double get_postprocessing_tolerance() { return this->postprocessing_tolerance; };
        float get_postprocessing_min_depth() { return this->postprocessing_min_depth; };
        double get_postprocessing_min_cosine_similarity() { return this->postprocessing_min_cosine_similarity; };
        size_t get_mst_min_cluster_size() { return this->mst_min_cluster_size; };
        bool get_no_postprocessing() { return this->no_postprocessing; };
        bool get_is4d() { return this->is4d; };
        double get_mst_max_mean_edge_len() { return this->mst_max_mean_edge_len; };
};

#endif