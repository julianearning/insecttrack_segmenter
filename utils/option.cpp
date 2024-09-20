#include "option.hpp"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>


Opt::Opt() {
    this->scale = 0.526315789;
    this->preprocessing_method = "octree";
    this->segmentation_method = "umst";
    this->smooth = false;
    this->smooth_radius = 35;
    this->octree_target_size = 5;
    this->octree_selection = "center";
    this->octree_denoise=1;
    this->eps=-1;
    this->minPts=6;
    this->max_lsq_distance=100;
    this->buffer_size=10;
    this->space_eps=40;
    this->min_cosine_score=0.999;
    this->mst_max_branch_depth=100;
    this->mst_max_edge_dist=-1;
    this->mst_twoD = false;
    this->mst_pca_angle = false;
    //this->mst_pca_angle_value = 0.005;
    this->mst_pca_angle_value = 0.005;
    this->mst_min_cluster_size=4;
    this->mst_max_mean_edge_len=-1;
    this->generate_preprocessing_debug=false;
    this->generate_preprocessing_debug_filename="debuggnuplot.gnuplot";
    this->out_filename = "test.csv";
    this->generate_mst_debug_gnuplot=false;
    this->mst_debug_gnuplot_filename = "debug_mst.gnuplot";
    this->remove_duplicates=false;
    this->window_mode=false;
    this->silent=false;
    this->postprocessing_threshold_sd=22.36;
    //this->postprocessing_min_depth=5;
    this->postprocessing_min_depth=0.04;
    this->postprocessing_point_proportion=0.3;
    this->postprocessing_tolerance=66.5*this->scale;
    this->postprocessing_min_cosine_similarity=0.7;
    this->no_postprocessing=false;
    this->is4d=false;
}


int Opt::parse_args(int argc, char ** argv ) {
    // parse command line
    std::pair<double, bool> tmp;
    try {
      for (int i = 2; i < argc; i++) {
        if (0 == strcmp(argv[i], "-scale")) {
            ++i;
            if (i < argc) {
                this->scale = std::stod(argv[i]);
            } else {
                return 1;
            }
        } else if (0 == strcmp(argv[i], "-octree_target_size")) {
            ++i;
            if (i < argc) {
                this->octree_target_size = std::stoi(argv[i]);
            } else {
                return 1;
            }
        } else if (0 == strcmp(argv[i], "-dbscan_eps")) {
            ++i;
            if (i < argc) {
                this->eps = std::stod(argv[i]);
            } else {
                return 1;
            }
        } else if (0 == strcmp(argv[i], "-postprocessing_threshold_sd")) {
            ++i;
            if (i < argc) {
                this->postprocessing_threshold_sd = std::stod(argv[i]);
            } else {
                return 1;
            }
        } else if (0 == strcmp(argv[i], "-postprocessing_min_depth")) {
            ++i;
            if (i < argc) {
                this->postprocessing_min_depth = (float)std::stod(argv[i]);
            } else {
                return 1;
            }
        } else if (0 == strcmp(argv[i], "-postprocessing_point_proportion")) {
            ++i;
            if (i < argc) {
                this->postprocessing_point_proportion = std::stod(argv[i]);
            } else {
                return 1;
            }
        } else if (0 == strcmp(argv[i], "-postprocessing_tolerance")) {
            ++i;
            if (i < argc) {
                this->postprocessing_tolerance = std::stod(argv[i]);
            } else {
                return 1;
            }
        } else if (0 == strcmp(argv[i], "-postprocessing_min_cosine_similarity")) {
            ++i;
            if (i < argc) {
                this->postprocessing_min_cosine_similarity = std::stod(argv[i]);
            } else {
                return 1;
            }
        } else if (0 == strcmp(argv[i], "-space_eps")) {
            ++i;
            if (i < argc) {
                this->space_eps = std::stod(argv[i]);
            } else {
                return 1;
            }
        } else if (0 == strcmp(argv[i], "-mst_max_edge_dist")) {
            ++i;
            if (i < argc) {
                this->mst_max_edge_dist = std::stod(argv[i]);
            } else {
                return 1;
            }
        } else if (0 == strcmp(argv[i], "-mst_max_branch_depth")) {
            ++i;
            if (i < argc) {
                this->mst_max_branch_depth = std::stoi(argv[i]);
            } else {
                return 1;
            }
        } else if (0 == strcmp(argv[i], "-mst_min_cluster_size")) {
            ++i;
            if (i < argc) {
                this->mst_min_cluster_size = std::stoi(argv[i]);
            } else {
                return 1;
            }
        } else if (0 == strcmp(argv[i], "-mst_max_mean_edge_len")) {
            ++i;
            if (i < argc) {
                this->mst_max_mean_edge_len = std::stod(argv[i]);
            } else {
                return 1;
            }
        } else if (0 == strcmp(argv[i], "-octree_denoise")) {
            ++i;
            if (i < argc) {
                this->octree_denoise = std::stoi(argv[i]);
            } else {
                return 1;
            }
        } else if (0 == strcmp(argv[i], "-out_filename")) {
            ++i;
            if (i < argc) {
                this->out_filename = argv[i];
            } else {
                return 1;
            }
        } else if (0 == strcmp(argv[i], "-dbscan_minPts")) {
            ++i;
            if (i < argc) {
                this->minPts = std::stoi(argv[i]);
            } else {
                return 1;
            }
        } else if (0 == strcmp(argv[i], "-buffer_size")) {
            ++i;
            if (i < argc) {
                this->buffer_size = std::stoi(argv[i]);
            } else {
                return 1;
            }
        } else if (0 == strcmp(argv[i], "-max_lsq_distance")) {
            ++i;
            if (i < argc) {
                this->max_lsq_distance = std::stod(argv[i]);
            } else {
                return 1;
            }
        } else if (0 == strcmp(argv[i], "-min_cosine_score")) {
            ++i;
            if (i < argc) {
                this->min_cosine_score = std::stod(argv[i]);
            } else {
                return 1;
            }
        } else if (0 == strcmp(argv[i], "-generate_preprocessing_debug")) {
            ++i;
            if (i < argc) {
                this->generate_preprocessing_debug = true;
                this->generate_preprocessing_debug_filename = std::string(argv[i]);
            } else {
                return 1;
            }
        } else if (0 == strcmp(argv[i], "-window_mode")) {
            this->window_mode = true;
        } else if (0 == strcmp(argv[i], "-silent")) {
            this->silent = true;
        } else if (0 == strcmp(argv[i], "-mst_twoD")) {
            this->mst_twoD = true;
        } else if (0 == strcmp(argv[i], "-4d")) {
            this->is4d = true;
        } else if (0 == strcmp(argv[i], "-mst_pca_angle")) {
            ++i;
            if (i < argc) {
                this->mst_pca_angle=true;
                this->mst_pca_angle_value=std::stod(argv[i]);
            } else {
                return 1;
            }
        } else if (0 == strcmp(argv[i], "-preprocessing_method")) {
          ++i;
          if (i < argc) {
            if((strcmp(argv[i], "octree")!=0) && (strcmp(argv[i], "space")!=0) && (strcmp(argv[i], "none")!=0)) {
              return 1;
            }
            this->preprocessing_method = std::string(argv[i]);
          } else {
            return 1;
          }
        } else if (0 == strcmp(argv[i], "-smooth")) {
            this->smooth = true;
        } else if (0 == strcmp(argv[i], "-no_postprocessing")) {
            this->no_postprocessing = true;
        } else if (0 == strcmp(argv[i], "-remove_duplicates")) {
            this->remove_duplicates = true;
        } else if (0 == strcmp(argv[i], "-smooth_radius")) {
            ++i;
            if (i < argc) {
                this->smooth_radius = std::stod(argv[i]);
            } else {
                return 1;
            }
        } else if (0 == strcmp(argv[i], "-generate_mst_debug_gnuplot")) {
            ++i;
            if (i < argc) {
                this->generate_mst_debug_gnuplot=true;
                this->mst_debug_gnuplot_filename=std::string(argv[i]);
            } else {
                return 1;
            }
        } else if (0 == strcmp(argv[i], "-segmentation_method")) {
          ++i;
          if (i < argc) {
            if((strcmp(argv[i], "dbscan")!=0) && (strcmp(argv[i], "umst")!=0) && (strcmp(argv[i], "dmst")!=0) && (strcmp(argv[i], "local_pca")!=0) && (strcmp(argv[i], "triplclust")!=0)) {
              return 1;
            }
            this->segmentation_method = std::string(argv[i]);
          } else {
            return 1;
          }
        } else if (0 == strcmp(argv[i], "-octree_selection")) {
          ++i;
          if (i < argc) {
            if((strcmp(argv[i], "centroid")!=0) && (strcmp(argv[i], "center")!=0)) {
              return 1;
            }
            this->octree_selection = std::string(argv[i]);
          } else {
            return 1;
          }
        } else if(0 == strcmp(argv[i], "-help")) { 
            return 1;
        } else {
          std::cout<<argv[i]<<": unknown argument!"<<std::endl;
          return 1;
        } 
      }
    } catch (const std::invalid_argument &e) {
      std::cerr << e.what() << std::endl;
      return 1;
    }
    return 0;
}
