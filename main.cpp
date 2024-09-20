#include <iostream>
#include <cmath>
#include <cstring>
#include <string>
#include <chrono>
#include <algorithm>
#include <cstdlib>
#include "utils/option.hpp"
#include "utils/fileio.hpp"
#include "utils/pointcloud.hpp"
#include "preprocessing/octree.hpp"
//#include "preprocessing/space.hpp"
#include "dbscan/dbscan.hpp"
//#include "kdtree/kdtree.hpp"
#include "mst/mst.hpp"
#include "assign_to_cluster/assign.hpp"
#include "loess/loess.hpp"
#include "postprocess/postprocess.hpp"
#include <flann/flann.hpp>


bool compareByT_temp(const Point& p1, const Point& p2) { 
    return p1.t < p2.t; 
}


//-------------------------------------------------------------------
// Smoothing of the PointCloud *cloud*.
// For every point the nearest neighbours in the radius *r* is searched
// and the centroid of this neighbours is computed. The result is
// returned in *result_cloud* and contains these centroids. The
// centroids are duplicated in the result cloud, so it has the same
// size and order as *cloud*.
//-------------------------------------------------------------------
void smoothen_cloud(PointCloud &cloud, double r) {
  PointCloud result_cloud;
  size_t n = cloud.size();

  // If the smooth-radius is zero return the unsmoothed pointcloud
  if (r == 0) {
    return;
  }

  size_t ndims=cloud[0].data.size()+1;
   
  float * data = new float[n*ndims];
  int temp;
  for(int i=0; i<n; i++) {
    temp=i*ndims;
    for(size_t j = 0; j<(ndims-1); j++) {
        *((data + temp)+j) = (float)cloud.at(i).data.at(j);  
    }
    *((data + temp)+(ndims-1)) = (float)cloud.at(i).t;
  }    
  flann::Matrix<float> dataset=flann::Matrix<float>(data, n, ndims);   // construct dataset 
  flann::Matrix<float> query=flann::Matrix<float>(data, n, ndims);     // construct queries (all datapoints again)     
  flann::SearchParams sp = flann::SearchParams(124);
  //flann::SearchParams sp = flann::SearchParams(32);
  sp.cores=0;  
  flann::Index<flann::L2<float> > index(dataset, flann::KDTreeIndexParams(32));  
  index.buildIndex();                                              // build KdTree with Index      
  // initialize datastructures for radiusSearch
  std::vector<std::vector<size_t>> indices;
  std::vector<std::vector<float>> dists;   
  // radius search for all points in pc
  index.radiusSearch(query, indices, dists, (float)r, sp);


  size_t result_size;
  Point new_point=Point(ndims);
  for(size_t i = 0; i<n; i++) {
    new_point=Point(ndims);
    result_size = indices.at(i).size();
    for(size_t j = 0; j<result_size; j++) {
        //std::cout<<new_point<<"   +   "<< cloud.at(indices.at(i).at(j));
        new_point = new_point + cloud.at(indices.at(i).at(j));
        //std::cout<<"   =   "<<new_point<<"\n";
    }
    new_point = new_point / result_size;
    result_cloud.push_back(new_point);
  }
  cloud.clear();
  for(size_t i = 0; i<result_cloud.size(); i++) {
    cloud.push_back(result_cloud.at(i));
  }
  std::sort(cloud.begin(), cloud.end(), compareByT_temp);
}

void remove_duplicates(PointCloud & cloud) {
    //if(!simple) std::sort(cloud.begin(), cloud.end(), [](Point a, Point b) { return a.t > b.t;});
    PointCloud output;
    for(size_t i = 1; i<cloud.size(); i++) {
        if(!(cloud.at(i-1)==cloud.at(i))) {
            output.push_back(cloud.at(i));
        } 
    }
    cloud.clear();
    for(size_t i = 0; i<output.size(); i++) {
        cloud.push_back(output.at(i));
    }
}


// usage message
const char *usage =
    "Usage:\n"
    "\tout <infile> [options]\n"
    "Options (defaults in brackets):\n"
    "\t-scale <scale>    Multiply t by scale. Default: 0.526315789=1/1.9\n"
    "\t-out_filename <filename> Filename for output csv\n"
    "\t-preprocessing_method [octree,space,none]  Preprocessing method\n"
    "\t-space_eps <eps> distance to consider for preprocessing method space\n"
    "\t-octree_target_size <size>   OcTree Target Size\n"
    "\t-octree_selection [center, centroid]  Choose point closest to what\n"    
    "\t-octree_denoise <number_of_points> Points that are in a box with less than or equal number_of_points points\n"
    "\t-smooth Smooth before segmenting\n" 
    "\t-smooth_radius <radius> Smoothing radius\n" 
    "\t-remove_duplicates Remove duplicates before preprocessing (good for MST?)\n"
    "\t-segmentation_method [dbscan,local_pca,umst,dmst]\n"
    "\t-dbscan_eps <eps> Eps parameter for dbscan. Choose -1 and preprocessing method octree for automatic choice. \n"
    "\t-dbscan_minPts <minPts> minPts. Default: 6. Reduce if tracks are seperated too much.\n"
    "\t-mst_max_edge_dist <max_edge_dist|-1> Max edge dist for MST method, choose -1 for automatic choice, choose -2 for inconsistent edge method\n"
    "\t-mst_max_branch_depth <max_branch_dist> Max branch dist for MST method\n"
    "\t-mst_pca_angle <cos_similarity_threshold> Use PCA angle to remove edges. Default: 0.005\n" 
    "\t-mst_min_cluster_size <min_cluster_size> All cluster with size<min_cluster_size are considered noise (label -1). Default: 2\n"
    "\t-mst_max_mean_edge_len <max_mean_edge_len> Helps with noise reduction with umst method. Default -1 (not used)\n"
    "\t-generate_mst_debug_gnuplot <filename> Depict the generated mst\n"
    "\t-max_lsq_distance <distance> for postprocessing??\n"
    "\t-min_cosine_score <score_threshold> for local pca method\n"
    "\t-buffer_size <buffer_size> for local pca method\n"
    "\t-generate_preprocessing_debug <filename>\n"
    "\t-window_mode Merge file with another file (overlapping)\n"
    "\t-silent Silent mode\n"
    "\t-4d\n"
    "\t-postprocessing_threshold_sd <threshold> The threshold for the standard deviation to classify compound trajectories (default: 500.0)\n"
    "\t-postprocessing_min_depth <min_depth> Minimum depth to detect a branch in the MST (default: 5)\n"
    "\t-postprocessing_point_proportion <point_proportion> the proportion of points that are used as begin and end segment (default: 0.3)\n"
    "\t-postprocessing_tolerance <tolerance> The tolerance in the t variable for segments to be candidates to merge (default: 100.0)\n"
    "\t-no_postprocessing Turn off postprocessing\n"
    "\t-postprocessing_min_cosine_similarity <min_cosine_similarity> Minumum cosine similarity (default: 0.5)\n"; 

int main(int argc, char * argv []) {


    std::vector<int> labels;
    std::vector<int> full_labels;
    PointCloud downsampled = PointCloud();
    std::vector<size_t> noise_points;


    std::chrono::steady_clock::time_point begin;
    std::chrono::steady_clock::time_point end;
    if(argc<2) {
        std::cout<<"you forgot the filename!"<<std::endl;
        return 1;
    }


    begin = std::chrono::steady_clock::now();
    Opt opt_params;
    if (opt_params.parse_args(argc, argv) != 0) {
        std::cerr << usage << std::endl;
        return 1;
    }
    end = std::chrono::steady_clock::now();
    bool silent=opt_params.get_silent();
    if(!silent) std::cout<<"Parsing Command Line Arguments: "<<std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()<<"ms\n\n";



    const size_t ndims = opt_params.get_is4d() ? 4 : 3;

    begin = std::chrono::steady_clock::now();

    PointCloud pc = PointCloud();
    // load csv file into pointcloud
    FileParser fp = FileParser();
    if(!silent) std::cout<<"Reading in file "<<std::string(argv[1])<<" with scale: "<<opt_params.get_scale()<<"\n";
    fp.parse_file(std::string(argv[1]), pc, opt_params.get_scale(), opt_params.get_is4d());

    end = std::chrono::steady_clock::now();
    if(!silent) std::cout<<"Reading data: "<<std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()<<"ms\n\n";

    ////// NO PREPROCESSING
    if(opt_params.get_preprocessing_method()=="none") {
        if(!silent) std::cout<<"No preprocessing\n";
        begin = std::chrono::steady_clock::now();
        std::vector<int> labels;
        size_t max_label=0;
        if(opt_params.get_segmentation_method() == "dbscan") {
            if(!silent) std::cout<<"Using segmentation method dbscan with eps: "<<opt_params.get_dbscan_eps()<<", and minPts: "<<opt_params.get_dbscan_minPts()<<"\n";
            // call DBSCAN on pointcloud
            DBSCAN db = DBSCAN(&pc);
            double eps = opt_params.get_dbscan_eps();
            if(opt_params.get_dbscan_eps()<0) {
                eps=22;
            }
            // labels eps minPts 
            db.calc_DBSCAN(&labels, eps*eps, opt_params.get_dbscan_minPts(), true, false, opt_params.get_buffer_size(), opt_params.get_max_lsq_distance()*opt_params.get_max_lsq_distance());
        } else if(opt_params.get_segmentation_method()== "local_pca") {
            if(!silent) std::cout<<"Using segmentation method local_pca with eps: "<<opt_params.get_dbscan_eps()<<", and minPts: "<<opt_params.get_dbscan_minPts()<<", and buffer_size: "<<opt_params.get_buffer_size()<<", and min_cosine_score: "<<opt_params.get_min_cosine_score()<<"\n";
            DBSCAN db = DBSCAN(&pc);
            db.calc_DBSCAN(&labels, opt_params.get_dbscan_eps()*opt_params.get_dbscan_eps(), opt_params.get_dbscan_minPts(), false, true, opt_params.get_buffer_size(), opt_params.get_min_cosine_score());
        } else if(opt_params.get_segmentation_method()=="umst") { 
            double eps_local=opt_params.get_mst_max_edge_dist();
            if((opt_params.get_mst_max_edge_dist() == -1)) {
                eps_local =50.0;
            } 
            if(!silent) std::cout<<"Using segmentation method undirected mst with max_edge_dist: "<<eps_local<<", and max_branch_depth: "<<opt_params.get_mst_max_branch_depth()<<"\n";
            UndirectedMST mst = UndirectedMST(&pc, &labels, (float)(eps_local*eps_local), opt_params.get_mst_min_cluster_size(), ndims, (opt_params.get_mst_max_mean_edge_len()>0 ? (opt_params.get_mst_max_mean_edge_len()*opt_params.get_mst_max_mean_edge_len()) : -1.0));
            max_label = mst.colorMST()-1;
            //const char * fname = "debug_mst.gnuplot";
            if(opt_params.get_generate_mst_debug_gnuplot()) {
               mst.draw_mst(opt_params.get_mst_debug_gnuplot_filename().c_str(), opt_params.get_is4d()); 
            }
        }  else if(opt_params.get_segmentation_method()=="dmst") {
            double eps_local=opt_params.get_mst_max_edge_dist();
            if((opt_params.get_mst_max_edge_dist() == -1)) {
                eps_local = 70.0;
            } 
            if(!silent) std::cout<<"Using segmentation method directed mst with max_edge_dist: "<<eps_local<<", and max_branch_depth: "<<opt_params.get_mst_max_branch_depth()<<"\n";
            DirectedMST mst = DirectedMST(&pc, &labels, (float)(eps_local*eps_local), opt_params.get_mst_min_cluster_size(), ndims);
            max_label = mst.colorMST()-1;
            //const char * fname = "debug_mst.gnuplot";
            if(opt_params.get_generate_mst_debug_gnuplot()) {
               mst.draw_mst(opt_params.get_mst_debug_gnuplot_filename().c_str(), opt_params.get_is4d()); 
            }
        }
        end = std::chrono::steady_clock::now();
        if(!silent) std::cout<<"Segmentation: "<<std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()<<"ms\n\n";



        /////////////// Postprocessing ///////////////////
        if(!opt_params.get_no_postprocessing()) {
            begin = std::chrono::steady_clock::now();
            PostProcessor pp = PostProcessor(&pc, &labels, (float)opt_params.get_postprocessing_threshold_sd(), opt_params.get_postprocessing_min_depth(), (size_t)max_label, (float)opt_params.get_postprocessing_point_proportion(), opt_params.get_postprocessing_tolerance(), opt_params.get_postprocessing_min_cosine_similarity(), opt_params.get_mst_min_cluster_size(), ndims);
            size_t num_labels_splitted = pp.postprocess();
            end = std::chrono::steady_clock::now();
            if(!silent) std::cout<<"Postprocessing: "<<std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()<<"ms\n\n";
        }


        begin = std::chrono::steady_clock::now();
        std::vector<int> label_order;
        std::vector<PointCloud> loess_vectors;
        LOESS loess = LOESS(&pc, &labels, &loess_vectors, &label_order, 0.2, ndims);
        end = std::chrono::steady_clock::now();
        if(!silent) std::cout<<"LOESS: "<<std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()<<"ms\n\n";


        begin = std::chrono::steady_clock::now();
        if(!silent) std::cout<<"Writing the pointcloud + labels"<<std::endl;
        if(opt_params.get_window_mode()) {
            fp.merge_csv(opt_params.get_out_filename(), pc, labels, opt_params.get_scale(), opt_params.get_is4d());
        } else {
            fp.write_csv(opt_params.get_out_filename(), pc, labels, false); 
        }
        end = std::chrono::steady_clock::now();
        if(!silent) std::cout<<"Writing: "<<std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()<<"ms\n\n";
        return 0;
    }


    ////////////////// PREPROCESSING //////////////////
    begin = std::chrono::steady_clock::now();

    if(opt_params.get_preprocessing_method() == "octree") {
        if(!silent) std::cout<<"Using preprocessing method octree with octree_target_size: "<<opt_params.get_octree_target_size()<<"\n";
        OcTree oc = OcTree(&pc, &downsampled, opt_params.get_octree_target_size(), opt_params.get_octree_denoise(), ndims); 
        oc.preprocess(&noise_points);
        if((!silent) && (opt_params.get_octree_denoise()>0)) std::cout<<"Removed "<<noise_points.size()<<" noise points\n";
        //if(!silent) std::cout<<"Size of cloud without preprocessing: "<<pc.size()<<", after preprocessing: "<<downsampled.size()<<"\n";
    } /* else if(opt_params.get_preprocessing_method() == "space") {
        if(!silent) std::cout<<"Using preprocessing method space with space_eps: "<<opt_params.get_space_eps()<<"\n";
        Space sp = Space(&pc, &downsampled, opt_params.get_space_eps()*opt_params.get_space_eps());   // uses flann library, so all distances are internally not squared
        sp.preprocess();
    } */

    //fp.write_csv_no_label("preprocessing.csv", downsampled);

    if(opt_params.get_smooth()) {
        //if(!silent) std::cout<<"Smoothening point cloud with smooth_radius: "<<opt_params.get_smooth_radius()<<"\n";
        //smoothen_cloud(downsampled, 10);
        //if(!silent) std::cout<<"Pointcloud smoothed"<<std::endl;

        if(!silent) std::cout<<"Smoothening point cloud with smooth_radius: "<<opt_params.get_smooth_radius()<<"\n";
        smoothen_cloud(downsampled, opt_params.get_smooth_radius()*opt_params.get_smooth_radius());
        if(!silent) std::cout<<"Pointcloud smoothed"<<std::endl;

        // TEST TODODDDODODODODOD

        //PointCloud downsampled2=downsampled;
        //downsampled.clear();
//
        //OcTree oc = OcTree(&downsampled2, &downsampled, opt_params.get_octree_target_size(),0);  
//
        //oc.preprocess();

    }

    //fp.write_csv_no_label("preprocessing_after_smooth_eps"+std::to_string(opt_params.get_space_eps())+".csv", downsampled);


    
    if(opt_params.get_preprocessing_method()!="none" && opt_params.get_remove_duplicates()) {
        if(!silent) std::cout<<"Removing duplicates...\n";
        remove_duplicates(downsampled);
        //if(!silent) std::cout<<"[Done]\n";
    }
    
    if(opt_params.get_preprocessing_method()!="none") {
        if(!silent) std::cout<<"Size of cloud without preprocessing: "<<pc.size()<<", after preprocessing: "<<downsampled.size()<<"\n";
        //fp.write_csv_no_label("preprocessing_smooth_eps40_bigboyloop.csv", downsampled);
    }

    end = std::chrono::steady_clock::now();
    if(!silent) std::cout<<"Preprocessing: "<<std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()<<"ms\n\n";

    

    if(opt_params.get_preprocessing_method()!="none" && opt_params.get_generate_preprocessing_debug()) {
        std::string filename_debug_cpp = opt_params.get_generate_preprocessing_debug_filename();
        if(!silent) std::cout<<"Generating preprocessing debug gnuplot file and writing to: "<<filename_debug_cpp<<"\n";
        const char * filename_debug = filename_debug_cpp.c_str();
        fp.debug_gnuplot(pc, downsampled, filename_debug, opt_params.get_is4d());
    }

    
    /////////////////// SEGMENTATION ///////////////////
    begin = std::chrono::steady_clock::now();
    int max_label=0;
    if(opt_params.get_segmentation_method() == "dbscan") {
        if(!silent) std::cout<<"Using segmentation method dbscan with eps: "<<opt_params.get_dbscan_eps()<<", and minPts: "<<opt_params.get_dbscan_minPts()<<"\n";
        // call DBSCAN on pointcloud
        DBSCAN db = DBSCAN(&downsampled);
        double eps_local=opt_params.get_dbscan_eps();
        if((opt_params.get_preprocessing_method() == "octree") && (opt_params.get_dbscan_eps()<0)) {
            //eps_local = 40;
            eps_local = 6*std::sqrt(2*opt_params.get_octree_target_size()*opt_params.get_octree_target_size());
            //eps_local = 8 * (double)opt_params.get_octree_target_size();
            if(!silent) std::cout<<"Automatic choice of dbscan_eps="<<eps_local<<"\n";
        } else if(opt_params.get_dbscan_eps()<0){
            if(!silent) std::cout<<"[ Warning ] Choosing dbscan_eps < 0 without using preprocessing method octree is impossible. Choosing fallback eps value: 40\n";
            eps_local = 40;
        }
        // labels eps minPts 
        max_label = db.calc_DBSCAN(&labels, eps_local*eps_local, opt_params.get_dbscan_minPts(), true, false, opt_params.get_buffer_size(), opt_params.get_max_lsq_distance());
    } else if(opt_params.get_segmentation_method()== "local_pca") {
        if(!silent) std::cout<<"Using segmentation method local_pca with eps: "<<opt_params.get_dbscan_eps()<<", and minPts: "<<opt_params.get_dbscan_minPts()<<", and buffer_size: "<<opt_params.get_buffer_size()<<", and min_cosine_score: "<<opt_params.get_min_cosine_score()<<"\n";
        DBSCAN db = DBSCAN(&downsampled);
        double eps_local=opt_params.get_dbscan_eps();
        if((opt_params.get_preprocessing_method() == "octree") && (opt_params.get_dbscan_eps()<0)) {
            //eps_local = 40;
            eps_local = 6*std::sqrt(2*opt_params.get_octree_target_size()*opt_params.get_octree_target_size());
            //eps_local = 8 * opt_params.get_octree_target_size();
            if(!silent) std::cout<<"Automatic choice of dbscan_eps="<<eps_local<<"\n";
        } else if(opt_params.get_dbscan_eps()<0){
            if(!silent) std::cout<<"[ Warning ] Choosing dbscan_eps < 0 without using preprocessing method octree is impossible. Choosing fallback eps value: 40\n";
            eps_local = 40;
        }
        max_label = db.calc_DBSCAN(&labels, eps_local*eps_local, opt_params.get_dbscan_minPts(), false, true, opt_params.get_buffer_size(), opt_params.get_min_cosine_score());
    } else if(opt_params.get_segmentation_method()=="umst") {
        //DirectedMST mst = DirectedMST(&downsampled, &labels, opt_params.get_mst_max_edge_dist()*opt_params.get_mst_max_edge_dist(), opt_params.get_mst_max_branch_depth(), true, opt_params.get_mst_pca_angle(), opt_params.get_mst_pca_angle_value(), opt_params.get_preprocessing_method() == "octree" ? opt_params.get_octree_target_size() : -1);
        double eps_local=opt_params.get_mst_max_edge_dist();
        if((opt_params.get_mst_max_edge_dist() == -1)) {
            if(opt_params.get_preprocessing_method() == "octree")
                eps_local = 5*std::sqrt(2*opt_params.get_octree_target_size()*opt_params.get_octree_target_size());
            else
                eps_local = 35.0;
        } 
        if(!silent) std::cout<<"Using segmentation method undirected mst with max_edge_dist: "<<eps_local<<", and max_branch_depth: "<<opt_params.get_mst_max_branch_depth()<<"\n";
        UndirectedMST mst = UndirectedMST(&downsampled, &labels, (float)(eps_local*eps_local), opt_params.get_mst_min_cluster_size(), ndims, (opt_params.get_mst_max_mean_edge_len()>0 ? (opt_params.get_mst_max_mean_edge_len()*opt_params.get_mst_max_mean_edge_len()) : -1.0));
        max_label = mst.colorMST()-1;
        //UndirectedMST mst = UndirectedMST(&downsampled, &labels, opt_params.get_mst_max_edge_dist()*opt_params.get_mst_max_edge_dist(), opt_params.get_mst_max_branch_depth(), true, opt_params.get_mst_pca_angle(), opt_params.get_mst_pca_angle_value(), opt_params.get_preprocessing_method() == "octree" ? opt_params.get_octree_target_size() : -1);
        //const char * fname = "debug_mst.gnuplot";
        if(opt_params.get_generate_mst_debug_gnuplot()) {
           mst.draw_mst(opt_params.get_mst_debug_gnuplot_filename().c_str(), opt_params.get_is4d());
        }
    } else if(opt_params.get_segmentation_method()=="dmst") {
        //DirectedMST mst = DirectedMST(&downsampled, &labels, opt_params.get_mst_max_edge_dist()*opt_params.get_mst_max_edge_dist(), opt_params.get_mst_max_branch_depth(), true, opt_params.get_mst_pca_angle(), opt_params.get_mst_pca_angle_value(), opt_params.get_preprocessing_method() == "octree" ? opt_params.get_octree_target_size() : -1);
        double eps_local=opt_params.get_mst_max_edge_dist();
        if((opt_params.get_mst_max_edge_dist() == -1)) {
            if(opt_params.get_preprocessing_method() == "octree")
                eps_local = 10*std::sqrt(2*opt_params.get_octree_target_size()*opt_params.get_octree_target_size());
            else
                eps_local = 70;
        } 
        if(!silent) std::cout<<"Using segmentation method directed mst with max_edge_dist: "<<eps_local<<", and max_branch_depth: "<<opt_params.get_mst_max_branch_depth()<<"\n";
        DirectedMST mst = DirectedMST(&downsampled, &labels, (float)(eps_local*eps_local), opt_params.get_mst_min_cluster_size(), ndims);
        max_label = mst.colorMST()-1;
        //UndirectedMST mst = UndirectedMST(&downsampled, &labels, opt_params.get_mst_max_edge_dist()*opt_params.get_mst_max_edge_dist(), opt_params.get_mst_max_branch_depth(), true, opt_params.get_mst_pca_angle(), opt_params.get_mst_pca_angle_value(), opt_params.get_preprocessing_method() == "octree" ? opt_params.get_octree_target_size() : -1);
        //const char * fname = "debug_mst.gnuplot";
        if(opt_params.get_generate_mst_debug_gnuplot()) {
           mst.draw_mst(opt_params.get_mst_debug_gnuplot_filename().c_str(), opt_params.get_is4d());
        }
    }
    end = std::chrono::steady_clock::now();
    if(!silent) std::cout<<"Segmentation: "<<std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()<<"ms\n\n";

    /////////////// Postprocessing ///////////////////
    if(!opt_params.get_no_postprocessing()) {
        begin = std::chrono::steady_clock::now();
        PostProcessor pp = PostProcessor(&downsampled, &labels, (float)opt_params.get_postprocessing_threshold_sd(), opt_params.get_postprocessing_min_depth(), (size_t)max_label, (float)opt_params.get_postprocessing_point_proportion(), opt_params.get_postprocessing_tolerance(), opt_params.get_postprocessing_min_cosine_similarity(), opt_params.get_mst_min_cluster_size(), ndims);
        size_t num_labels_splitted = pp.postprocess();
        end = std::chrono::steady_clock::now();
        if(!silent) std::cout<<"Postprocessing: "<<std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()<<"ms\n\n";
    }

    
    //begin = std::chrono::steady_clock::now();
    //std::vector<PointCloud> loess_vectors;
    //std::vector<int> label_order;
    //LOESS loess = LOESS(&downsampled, &labels, &loess_vectors, &label_order, 3);
    //end = std::chrono::steady_clock::now();
    //if(!silent) std::cout<<"LOESS: "<<std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()<<"ms\n\n";


    /////////////////// Assigning ///////////////////
    begin = std::chrono::steady_clock::now();
    if(!silent) std::cout<<"Assigning unlabelled points to closest cluster"<<std::endl;
    if(opt_params.get_preprocessing_method() == "octree") {
        Assign(&downsampled, &pc, &labels, &full_labels, &noise_points, ndims);
    } else {
        Assign(&downsampled, &pc, &labels, &full_labels, ndims);
    }
    end = std::chrono::steady_clock::now();
    if(!silent) std::cout<<"Assigning: "<<std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()<<"ms\n\n";


    /////////////////// LOESS ///////////////////

    begin = std::chrono::steady_clock::now();
    std::vector<PointCloud> loess_vectors;
    std::vector<int> label_order;
    LOESS loess = LOESS(&downsampled, &labels, &loess_vectors, &label_order, 0.2, ndims);
    /* if(opt_params.get_smooth_radius()>0) {
        LOESS loess = LOESS(&downsampled, &labels, &loess_vectors, &label_order, 3,opt_params.get_smooth_radius()*opt_params.get_smooth_radius());
    } else {
        LOESS loess = LOESS(&downsampled, &labels, &loess_vectors, &label_order, 3,-1.0);
    } */
    end = std::chrono::steady_clock::now();
    if(!silent) std::cout<<"LOESS: "<<std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()<<"ms\n\n";


    /////////////////// Writing ///////////////////
    begin = std::chrono::steady_clock::now();
    if(!silent) std::cout<<"Writing the pointcloud + labels"<<std::endl;
    if(opt_params.get_window_mode()) {
        fp.merge_csv(opt_params.get_out_filename(), pc, full_labels, opt_params.get_scale(), opt_params.get_is4d());
    } else {
        fp.write_csv(opt_params.get_out_filename(), pc, full_labels, false, &noise_points); 
        fp.write_csv("downsampled.csv", downsampled, labels, false, &noise_points); 
        fp.write_loess("loess.csv", &loess_vectors, &label_order, false, opt_params.get_scale());
    }
    end = std::chrono::steady_clock::now(); 
    if(!silent) std::cout<<"Writing: "<<std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()<<"ms\n\n";   
    
    return 0;
}
