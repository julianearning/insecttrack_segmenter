#include "postprocess.hpp"

#include <cmath>
#include <iostream>
#include <unordered_map>
#include <boost/functional/hash.hpp>
#include "../mst/mst.hpp"
#include "../utils/mathstuff.hpp"
//#include "../utils/fileio.hpp"

void RANSAC(PointCloud * neighbourhood, float t, size_t * i, size_t * j, size_t ndims) {
    MathStuff ms;
    size_t max_i=0;
    size_t max_j=0;
    size_t max_num_supporters=0;
    size_t curr_supporters=0;
    Point a=Point(ndims);
    Point b=Point(ndims);
    Point c=Point(ndims);
    for(size_t i = 0; i<neighbourhood->size(); i++) {
        a=neighbourhood->at(i);
        for(size_t j = 0; j<i; j++) {
            b=neighbourhood->at(j)-neighbourhood->at(i);
            curr_supporters=0;
            for(size_t k = 0; k<neighbourhood->size(); k++) {
                c=neighbourhood->at(k);
                //std::cout<<ms.orthogonal_lsq_distance(&a,&b,&c)<<"\n";
                if(ms.orthogonal_lsq_distance(&a,&b,&c, ndims)<t) curr_supporters++;
            }
            if(curr_supporters>max_num_supporters) {
                max_i=i;
                max_j=j;
                max_num_supporters=curr_supporters;
            }
        }
        for(size_t j = i+1; j<neighbourhood->size(); j++) {
            b=neighbourhood->at(j)-neighbourhood->at(i);
            curr_supporters=0;
            for(size_t k = 0; k<neighbourhood->size(); k++) {
                c=neighbourhood->at(k);
                //std::cout<<ms.orthogonal_lsq_distance(&a,&b,&c)<<"\n";
                if(ms.orthogonal_lsq_distance(&a,&b,&c,ndims)<t) curr_supporters++;
            }
            if(curr_supporters>=max_num_supporters) {
                max_i=i;
                max_j=j;
                max_num_supporters=curr_supporters;
            }
        }
    }
    *i=max_i;
    *j=max_j;
}


void get_velocity(PointCloud * cloud, Point * start, Point * velocity) {
    size_t ndims = cloud->at(0).data.size()+1;
    Point a=Point(ndims);
    for(size_t i = 0; i<(int)(cloud->size()/2); i++) {
        a=a+Point(cloud->at(i));
    }
    a=a/(int)(cloud->size()/2);
    Point b=Point(ndims);
    for(size_t i = (int)(cloud->size()/2)+1; i<cloud->size(); i++) {
        b=b+Point(cloud->at(i));
    }
    b=b/(int)(cloud->size()/2);
    *start = b;
    *velocity = b-a;
}


double mymax(double a, double b) {
    return a>b?a:b;
}
double mymin(double a,double b) {
    return a<b?a:b;
}


double get_velocity_distance(Point * start, Point * velocity, Point * query) {
    // position an der query sein müsste bei seinem t wenn start und velocity gegeben ist, ist: 
    // p(t_q) = start + (t_q/(t_v-t_start) * v)
    double t_factor = query->t/(velocity->t-start->t);
    Point p_tq = (*start) + (*velocity) * t_factor;  // * for Point and scalar only works in one direction, sorry
    return query->spacial_euclidian_distance(p_tq);
}


PostProcessor::PostProcessor(PointCloud * pc, std::vector<int> * labels, float threshold_sd, float min_depth, size_t max_label, float point_proportion, double tolerance, double min_cosine_similarity, size_t min_cluster_size, size_t ndims) {
    this->pc = pc;
    this->labels = labels;
    this->threshold_sd = threshold_sd;
    this->min_depth = min_depth;
    this->max_label = max_label;
    this->point_proportion = point_proportion;
    this->tolerance = tolerance;
    this->min_cosine_similarity = min_cosine_similarity;
    this->min_cluster_size = min_cluster_size;
    this->ndims = ndims;
}


void PostProcessor::get_cluster(size_t c, PointCloud* input, std::vector<int>* input_labels, PointCloud* cluster) {
    cluster->clear();
    for(size_t i = 0; i<input->size(); i++) {
        if(input_labels->at(i)==c) {
            cluster->push_back(input->at(i)); 
        }
    }
}


float PostProcessor::compute_spacial_variance(PointCloud* cluster) {
    std::vector<float> diffs;
    Point a; 
    Point b;
    float mean=0;
    float sd=0;
    float dist;
    for(size_t i = 1; i<cluster->size(); i++) {
        a = cluster->at(i-1);
        b = cluster->at(i);
        diffs.push_back(a.spacial_euclidian_distance(b));
        mean+=dist;
    }
    for(size_t i = 0; i<diffs.size(); i++) {
        dist = diffs.at(i);
        sd+=(dist-mean)*(dist-mean);
    }
    // standard deviation ( except the square root is missing (maybe higher differences are weighted higher?) )
    sd = std::sqrt(sd/(cluster->size()-1));
    return sd;
}


/* void PostProcessor::get_velocity_of_cluster(std::pair<PointCloud,PointCloud>* segments, Point *a_end, Point *b_end) {
    get_velocity(&segments->second, b_end, a_end);
} */


void PostProcessor::get_pca_of_cluster(std::pair<PointCloud,PointCloud>* segments, Point *a_begin, Point *b_begin, Point *a_end, Point *b_end) {
    //MathStuff ms = MathStuff();
    size_t ii;
    size_t jj;
    //ms.orthogonal_lsq(segments->first, a_begin, b_begin);
    RANSAC(&segments->first, 100.0, &ii, &jj, ndims);
    //std::cout<<cluster->size()<<" "<<ii<<" "<<jj<<"\n";
    //*a_begin = segments->second.at(ii);
    *a_begin=segments->first.at(0);
    *b_begin = segments->second.at(jj)-segments->second.at(ii);

    //ms.orthogonal_lsq(segments->second, a_end, b_end);
    RANSAC(&segments->second, 100.0, &ii, &jj, ndims);
    //*a_end = segments->second.at(ii);
    *a_end=segments->second.at(segments->second.size()-1);
    *b_end = segments->second.at(jj)-segments->second.at(ii);
    //const size_t number_points = (size_t)(point_proportion * cluster->size());
    //*a_begin=sgmeents->first.at(0);
    //*b_begin=segments->first.at(0)-segments->first.at(segments->first.size()-1);
    //*a_end=segments->second.at(segments->second.size()-1);
    //*b_end=segments->second.at(segments->second.size()-1)-segments->second.at(0);
}


//
size_t PostProcessor::postprocess() {
    const size_t available_labels=max_label; // copy because max_label is getting changed
    PointCloud cluster;     // for points within one cluster
    float sd=0;
    size_t cool_label_index=0;
    size_t max=0;
    size_t num_seperated_points=0;
    size_t num_labels=0;
    Point a_begin=Point(ndims), b_begin=Point(ndims), a_end=Point(ndims), b_end=Point(ndims);
    std::unordered_map<size_t, size_t> n_extension_points; // for adding the subtrajectories together 
    std::vector<std::pair<PointCloud, PointCloud>> segments;
    //std::vector<std::pair<std::pair<Point,Point>, std::pair<Point, Point>>> pca_segments;
    PointCloud beginning;
    PointCloud end;
    PointCloud segment;
    size_t number_points = 10;
    size_t min_number_points=10;
    double distance=0.0;
    MathStuff ms = MathStuff();
    double time_difference=0.0;
    //double time_overlap=0.0;
    //std::vector<std::pair<std::pair<size_t, size_t>, std::pair<size_t, double>>> least_distances;
    //std::unordered_map<std::pair<size_t, size_t>,std::pair<size_t, double>, boost::hash<std::pair<size_t, size_t>>> least_distances;
    std::unordered_map<std::pair<size_t, size_t>, double, boost::hash<std::pair<size_t, size_t>>> least_distances;
    std::unordered_map<size_t, bool> begin_already_paired;
    std::pair<size_t, size_t> pair;
    std::vector<size_t> root_segments;
    //std::vector<size_t> labels_merged;   // stores 
    size_t max_i=0;
    //size_t min_k=0;
    size_t max_j=0;
    size_t curr_min_depth=this->min_depth;
    double max_distance=0.0;
    size_t curr_root_segment=0;

    for(size_t c=0; c<=available_labels; c++) { // iterate over all labels
        // collect the Points of label c
        get_cluster(c, pc, labels, &cluster);
        sd=0;

        // calculate the spacial distance between points
        sd = compute_spacial_variance(&cluster);
        //std::cout<<"sd: "<<sd<<" >= "<<threshold_sd<<"\n";

        // classifies as compound cluster?
        if(sd>=threshold_sd) {
            if(min_depth<1.0) {
                curr_min_depth = cluster.size()*this->min_depth;
                curr_min_depth = curr_min_depth > 30 ? 30 : curr_min_depth;
                //curr_min_depth = curr_min_depth >= cluster.size() ? cluster.size()*this->min_depth : curr_min_depth;
            }
            // seperate by detecting branches of sufficient length
            std::vector<int> temp_labels;
            UndirectedMST mst = UndirectedMST(&cluster, &temp_labels, curr_min_depth, min_cluster_size, ndims);
            num_labels = mst.colorMST();
            // branches detected? correct oversegmentation next
            if(num_labels > 2 ) {
                //FileParser fp = FileParser();
                //fp.write_csv("oversegmentation.csv", cluster, temp_labels, false);    
                segments.clear();
                least_distances.clear();
                begin_already_paired.clear();
                root_segments.clear();

                ///// Collect the beginning and end segments of each cluster (i.e. the clusters that were created by detecting branches)
                // store beginning and end segments of the clusters in segments datastructure
                for(size_t i = 0; i<num_labels; i++) {
                    root_segments.push_back(i);
                    get_cluster(i, &cluster, &temp_labels, &segment);
                    // number_points that are included in end and beginning are either a proportion of the number of points in the cluster, but never less than min_number_points and never more than the size of the cluster alltogether
                    number_points = (size_t)(point_proportion * segment.size());
                    number_points = number_points >= min_number_points ? number_points : min_number_points; 
                    number_points = number_points > segment.size() ? segment.size() : number_points;

                    segments.push_back(std::make_pair(PointCloud(),PointCloud()));
                    for(size_t j = 0; j<number_points; j++) {
                        segments.at(segments.size()-1).first.push_back(segment.at(j));
                    }
                    for(size_t j = (segment.size()-number_points); j<segment.size(); j++) {
                        segments.at(segments.size()-1).second.push_back(segment.at(j));
                    }
                }

                for(size_t i = 0; i<segments.size(); i++) {
                    begin_already_paired[i] = false;
                    for(size_t j = 0; j<segments.size(); j++) {
                        least_distances[std::make_pair(i,j)] = -1.0;  //std::make_pair(0, -1.0);
                    }
                }

                for(size_t i = 0; i<segments.size(); i++) {
                    end = segments.at(i).second;  // get end piece of segment i
                    ms.orthogonal_lsq(end, &a_end, &b_end, ndims); // calculate direction of beginning piece of segment j
                    for(size_t j = 0; j<segments.size(); j++) {
                        if(i!=j) {
                            //get_velocity(&segments.at(i).first, &a_begin, &b_begin);
                            time_difference = segments.at(j).first.at(0).t - end.at(end.size()-1).t;
                            if((time_difference>(-0.5*tolerance)) && (time_difference < tolerance)) {
                                ms.orthogonal_lsq(segments.at(j).first, &a_begin, &b_begin , ndims); // calculate direction of beginning piece of segment j

                                distance = b_end.cosine_similarity(b_begin);
                                pair = std::make_pair(i,j);
                                least_distances[pair]=std::abs(distance);

                            }
                        }
                    }
                }

                /* for(size_t i = 0; i<segments.size(); i++) {
                    for(size_t j = 0; j<segments.size(); j++) {
                        if(least_distances[std::make_pair(i,j)]!=-1) std::cout<<i<<", "<<j<<": "<<least_distances[std::make_pair(i,j)]<<"\n";  //std::make_pair(0, -1.0);
                    }
                } */

                for(size_t i = 0; i<(segments.size()-1); i++) {
                    max_distance=-1;
                    for(size_t j = 0; j<segments.size(); j++) {
                        if((least_distances[std::make_pair(i,j)]>max_distance) && (least_distances[std::make_pair(i,j)] >= 0)) {
                            max_distance = least_distances[std::make_pair(i,j)];
                            max_j=j;
                            max_i=i;
                        }
                    }
                    //std::cout<<i<<" "<<max_j<<" "<<max_distance<<"\n";
                    if(std::isfinite(max_distance) && (max_distance > min_cosine_similarity) && (!begin_already_paired[max_j])) {
                    //if(std::isfinite(max_distance) && (max_distance > min_cosine_similarity)) {
                        begin_already_paired[max_j] = true;
                        //std::cout<<"\t\tMerging cluster "<<max_i<<" and "<<max_j<<" with distance: "<<max_distance<<"\n";
                        // merge j and min_j, which means i have to merge root_segments(j) and root_segments(min_j) or the other way around when root_segments(j)>root_segments(min_j)
                        if(root_segments.at(max_i)<root_segments.at(max_j)) {
                            for(size_t k = 0; k<temp_labels.size(); k++) {
                                if(temp_labels.at(k)==root_segments.at(max_j)) temp_labels.at(k) = root_segments.at(max_i);
                            }
                            curr_root_segment=root_segments.at(max_j);
                            for(size_t k = 0; k<root_segments.size(); k++) {
                                if(root_segments.at(k)==root_segments.at(max_j)) root_segments.at(max_j) = root_segments.at(max_i);
                            }
                        } else {
                            for(size_t k = 0; k<temp_labels.size(); k++) {
                                if(temp_labels.at(k)==root_segments.at(max_i)) temp_labels.at(k) = root_segments.at(max_j);
                            }
                            curr_root_segment=root_segments.at(max_i);
                            for(size_t k = 0; k<root_segments.size(); k++) {
                                if(root_segments.at(k)==root_segments.at(max_i)) root_segments.at(max_i) = root_segments.at(max_j);
                            }
                        }
                    }
                }

                int max_a = 0;
                for(size_t k = 0; k<temp_labels.size(); k++) {
                    if(temp_labels.at(k)>max_a) max_a = temp_labels.at(k);
                }                


                // assign new labels
                cool_label_index=0;
                max=0;
                for(size_t i = 0; i<pc->size(); i++) {
                    if(labels->at(i)==c) {
                        if(temp_labels.at(cool_label_index)!=0) {
                            num_seperated_points++;
                            labels->at(i)=max_label+temp_labels.at(cool_label_index);
                            if(temp_labels.at(cool_label_index) > max) {
                                max=temp_labels.at(cool_label_index);
                            }
                        }
                        cool_label_index++;
                    }
                }
                max_label+=max;
                //mst.draw_mst((std::string("cluster_mst")+std::to_string(c)+std::string(".gnuplot")).c_str(), ndims==4);
            } 
        }
    }
    return num_seperated_points;
}