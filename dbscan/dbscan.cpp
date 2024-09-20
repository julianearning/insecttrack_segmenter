#include <flann/flann.hpp>
#include "dbscan.hpp"
#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include "../utils/mathstuff.hpp"


#define NOISE -1
#define UNDEFINED -2


DBSCAN::DBSCAN(PointCloud * pc_ptr) {
    this->pc_ptr = pc_ptr;
}

//
// orthogonal least squares fit with libeigen
//   pc: points
//   a, b: line representation as a + t*b
//   RC: largest eigenvalue
//
/* double DBSCAN::orthogonal_lsq(PointCloud & pc, Point* a, Point* b){
  double rc = 0.0;

  if (pc.size() == 0)
    return rc;
  
  // anchor point is mean value
  a->x = a->y = a->t = 0.0;
  for (size_t i=0; i < pc.size(); i++) {
    *a = *a + pc.at(i);
  }
  *a = *a / pc.size();
  
  // copy points to libeigen matrix
  Eigen::MatrixXf points = Eigen::MatrixXf::Constant(pc.size(), 3, 0);
  for (unsigned int i = 0; i < points.rows(); i++) {
    points(i,0) = pc.at(i).x;
    points(i,1) = pc.at(i).y;
    points(i,2) = pc.at(i).t;
  }

  // compute scatter matrix ...
  MatrixXf centered = points.rowwise() - points.colwise().mean();
  MatrixXf scatter = (centered.adjoint() * centered);

  // ... and its eigenvalues and eigenvectors
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> eig(scatter);
  Eigen::MatrixXf eigvecs = eig.eigenvectors();

  // we need eigenvector to largest eigenvalue
  // libeigen yields it as LAST column
  b->x = eigvecs(0,2); b->y = eigvecs(1,2); b->t = eigvecs(2,2);
  rc = eig.eigenvalues()(2);

  return (rc);
}
 

double orthogonal_lsq_distance(Point* a, Point* b, Point * c) {
  double distance = 0;
  double lambda;

  // copy points to libeigen vektor
  Eigen::Vector3f a_vec = Eigen::Vector3f(a->x,a->y,a->t);
  Eigen::Vector3f b_vec = Eigen::Vector3f(b->x,b->y,b->t);
  Eigen::Vector3f c_vec = Eigen::Vector3f(c->x,c->y,c->t);

  lambda = b_vec.transpose() * (c_vec - a_vec);

  distance = (c_vec - (a_vec + lambda * b_vec)).norm();

  //std::cout<<*a<<std::endl;
  //std::cout<<*b<<std::endl;
  //std::cout<<*c<<std::endl;
  //std::cout<<"------------------"<<std::endl;

  return distance; 
}


float angle(Point * a, Point * b) {
    return std::acos((float)(a->x*b->x+a->y*b->y+a->t*b->t)/(std::sqrt(a->x*a->x+a->y*a->y+a->t*a->t)*std::sqrt(b->x*b->x+b->y*b->y+b->t*b->t)));
}
*/


int DBSCAN::limited_calc_DBSCAN(std::vector<int> * labels_out, const double eps, const int minPts, const int buffer_size, const double max_lsq_distance) {

    //const int buffer_pts = 20;
    const bool t_limited = false;
    const bool angle_limited =false;

    const int n = pc_ptr->size();

    // initialize labels
    labels_out->clear();
    for(int i = 0; i<n; i++) {
        labels_out->push_back(UNDEFINED);
    }


    size_t ndims=pc_ptr->at(0).data.size()+1;
   
    float * data = new float[n*ndims];
    int temp;
    for(int i=0; i<n; i++) {
      temp=i*ndims;
      for(size_t j = 0; j<(ndims-1); j++) {
          *((data + temp)+j) = (float)pc_ptr->at(i).data.at(j);  
      }
      *((data + temp)+(ndims-1)) = (float)pc_ptr->at(i).t;
    }    
    flann::Matrix<float> dataset=flann::Matrix<float>(data, n, ndims);   // construct dataset 
    flann::Matrix<float> query=flann::Matrix<float>(data, n, ndims);     // construct queries (all datapoints again)     
    flann::SearchParams sp = flann::SearchParams(124);
    //flann::SearchParams sp = flann::SearchParams(32);
    sp.cores=0;  
    flann::Index<flann::L2<float> > index(dataset, flann::KDTreeIndexParams(4));  
    index.buildIndex();                                              // build KdTree with Index      
    // initialize datastructures for radiusSearch
    std::vector<std::vector<size_t>> indices;
    std::vector<std::vector<float>> dists;   
    // radius search for all points in pc
    index.radiusSearch(query, indices, dists, (float)eps, sp);



    // calculate how many neighbors have been returned
    int * n_neighbors = new int[n];
    for(int i = 0; i<n; i++) {
        n_neighbors[i] = indices.at(i).size();
    }
    
    // various stuff
    int status=0;
    // DBSCAN algorithm
    unsigned int c = 0;
    std::vector<int> s;
    // !!!!!!! adjusted part !!!!!!!!!!!
    int size_cluster = 0; 
    PointCloud buffer;
    buffer.reserve(buffer_size);
    double lsq_dist;
    Point buffer_pt=Point(ndims);
    Point b=Point(ndims);
    Point last_direction=Point(ndims);
    Point test1=Point(ndims);
    Point test2=Point(ndims);
    MathStuff mathstuff;
    //Point middle;
    //Point removed_point;
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // start DBSCAN algorithm
    for(int i = 0; i<n; ++i) {                          // iterate through all points (which are in chronological order)
        if(labels_out->at(i) == UNDEFINED) { 
            if(n_neighbors[i] >= minPts) {              // Core Point?
                c++;
                labels_out->at(i) = c;
                size_cluster = 1;
                buffer.clear();
                for(int j = 0; j<n_neighbors[i]; j++) { // add points to seedset s
                    s.push_back(indices.at(i).at(j));
                }
                while(! s.empty()) {                    // expand cluster
                    temp = s.back();
                    s.pop_back();

                    size_cluster++;
                    
                    // !!!!!!! adjusted part !!!!!!!!!!!
                    // update direction
                    if(buffer.size()>=buffer_size) {  // active phase
                        //removed_point = buffer.at(buffer.size()-1);
                        //removed_point/=buffer.size();
                        buffer.pop_back();
                        Point buffer_pt=Point(pc_ptr->at(temp));
                        //buffer_pt.x = pc_ptr->at(temp).x;
                        //buffer_pt.y = pc_ptr->at(temp).y;
                        //buffer_pt.t = pc_ptr->at(temp).t;
                        buffer.insert(buffer.begin(), buffer_pt);
                        mathstuff.orthogonal_lsq(buffer, &test1, &test2, ndims);
                        //middle = middle - removed_point + buffer_pt/buffer_size;
                    } else { 
                        // initialising phase
                        Point buffer_pt=Point(pc_ptr->at(temp));
                        //buffer_pt.x = pc_ptr->at(temp).x;
                        //buffer_pt.y = pc_ptr->at(temp).y;
                        //buffer_pt.t = pc_ptr->at(temp).t;
                        buffer.push_back(buffer_pt);
                        //middle += buffer_pt/buffer_size;
                    }

                    if(labels_out->at(temp) == NOISE) {
                        labels_out->at(temp) = c;
                    } 
                    if(labels_out->at(temp) == UNDEFINED) {
                        labels_out->at(temp) = c;
                        // !!!!!!! adjusted part !!!!!!!!!!!
                        if(buffer.size()>=buffer_size) {
                            Point c_pt=Point(pc_ptr->at(temp));
                            //c_pt.x = pc_ptr->at(temp).x;
                            //c_pt.y = pc_ptr->at(temp).y;
                            //c_pt.t = pc_ptr->at(temp).t;

                            //lsq_dist = mathstuff.orthogonal_lsq_distance(&test1, &test2, &c_pt);
                            lsq_dist = test2.cosine_similarity(test1-c_pt);

                            if((n_neighbors[temp]>=minPts) && (std::abs(lsq_dist) > max_lsq_distance)) {
                                for(int j = 0; j < n_neighbors[temp]; j++) {
                                    s.push_back(indices.at(temp).at(j));
                                }
                            }


                        } // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        else {
                            if(n_neighbors[temp]>=minPts) {
                                for(int j = 0; j < n_neighbors[temp]; j++) {
                                    //if((pc_ptr->at(i).t)<=pc_ptr->at(indices.at(temp).at(j)).t) {
                                        s.push_back(indices.at(temp).at(j));
                                    //}
                                }
                            }
                        }
                    }
                }
            } else {
                labels_out->at(i) = NOISE;
            }
        }
    }

    // cleaning up 
    delete [] n_neighbors;
    delete[] data;
    
    return c;
}




int DBSCAN::vanilla_calc_DBSCAN(std::vector<int> * labels_out, const double eps, const int minPts) {
    const int n = pc_ptr->size();

    // initialize labels
    labels_out->clear();
    for(int i = 0; i<n; i++) {
        labels_out->push_back(UNDEFINED);
    }
    // build indexed KDTREE here
    // first write data into helper array
    size_t ndims=pc_ptr->at(0).data.size()+1;
   
    float * data = new float[n*ndims];
    int temp;
    for(int i=0; i<n; i++) {
      temp=i*ndims;
      for(size_t j = 0; j<(ndims-1); j++) {
          *((data + temp)+j) = (float)pc_ptr->at(i).data.at(j);  
      }
      *((data + temp)+(ndims-1)) = (float)pc_ptr->at(i).t;
    }    
    flann::Matrix<float> dataset=flann::Matrix<float>(data, n, ndims);   // construct dataset 
    flann::Matrix<float> query=flann::Matrix<float>(data, n, ndims);     // construct queries (all datapoints again)     
    flann::SearchParams sp = flann::SearchParams(124);
    //flann::SearchParams sp = flann::SearchParams(32);
    sp.cores=0;  
    flann::Index<flann::L2<float> > index(dataset, flann::KDTreeIndexParams(4));  
    index.buildIndex();                                              // build KdTree with Index      
    // initialize datastructures for radiusSearch
    std::vector<std::vector<size_t>> indices;
    std::vector<std::vector<float>> dists;   
    // radius search for all points in pc
    index.radiusSearch(query, indices, dists, (float)eps, sp);

    /* double mean=0.0;
    int count=0;
    std::vector<float> test;
    for(auto i : dists) {
        if(i.size()>=2) {
            mean+=i.at(1);
            count++;
            test.push_back(i.at(1));
        }
    }
    std::cout<<mean/count<<std::endl;
    std::cout<<test.at((int)(test.size()/2))<<std::endl; */

    // calculate how many neighbors have been returned
    int * n_neighbors = new int[n];
    for(int i = 0; i<n; i++) {
        n_neighbors[i] = indices.at(i).size();
    }
    
    // DBSCAN algorithm
    unsigned int c = 0;      // cluster label
    int status=0;
    std::vector<int> s;
    int cnt=0;
    double last=pc_ptr->at(0).t;
    temp=0;
    for(int i = 0; i<n; ++i) {                          // iterate through all points (which are in chronological order)
        if(labels_out->at(i) == UNDEFINED) { 
            if(n_neighbors[i] >= minPts) {              // Core Point?
                c++;
                labels_out->at(i) = c;
                for(int j = 0; j<n_neighbors[i]; j++) { // add points to seedset s
                    s.push_back(indices.at(i).at(j));
                }
                while(! s.empty()) {                    // expand cluster
                    last=pc_ptr->at(temp).t;
                    temp = s.back();
                    s.pop_back();
                    if(labels_out->at(temp) == NOISE) {
                        labels_out->at(temp) = c;
                    } 
                    if(labels_out->at(temp) == UNDEFINED) {
                        labels_out->at(temp) = c;
                        if((n_neighbors[temp]>=minPts)) {
                                for(int j = 0; j < n_neighbors[temp]; j++) {
                                    s.push_back(indices.at(temp).at(j));
                                }
                        }
                    }
                }
            } else {
                labels_out->at(i) = NOISE;
            }
        }
    }

    // cleaning up
    delete [] n_neighbors;
    delete[] data;
    
    return c;
}



int DBSCAN::calc_DBSCAN(std::vector<int> * labels_out, const double eps, const int minPts, const bool vanilla, const bool limited, const int buffer_size, const double max_lsq_distance) {
    if(vanilla==limited) {
        std::cout<<"Either vanilla or limited flag must be true but not both"<<std::endl;
    }  
    if(vanilla) return vanilla_calc_DBSCAN(labels_out, eps, minPts);  
    if(limited) return limited_calc_DBSCAN(labels_out, eps, minPts, buffer_size, max_lsq_distance);
    return 0;
}





