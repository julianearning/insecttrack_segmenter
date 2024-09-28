#include "loess.hpp"
#include "WeightedLowess/WeightedLowess.hpp"
#include <iostream>
#include <algorithm>
#include <map>
#include <flann/flann.hpp>


extern bool compareByT(const Point0& p1, const Point0& p2);


//-------------------------------------------------------------------
// Smoothing of the PointCloud *cloud*.
// For every point the nearest neighbours in the radius *r* is searched
// and the centroid of this neighbours is computed. The result is
// returned in *result_cloud* and contains these centroids. The
// centroids are duplicated in the result cloud, so it has the same
// size and order as *cloud*.
//-------------------------------------------------------------------
void smoothen_cluster(PointCloud &cloud, double r) {
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
  Point0 new_point=Point0(ndims);
  for(size_t i = 0; i<n; i++) {
    new_point=Point0(ndims);
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
}


/* int find_and_return_index(std::vector<int> * labelled_points, int label) {
    for(size_t i = 0; i<labelled_points->size(); i++) {
        if(labelled_points->at(i) == label) return i;
    }
    return -1;
} */


LOESS::LOESS(PointCloud * pc, std::vector<int> * labels, std::vector<PointCloud> * loess, std::vector<int> * label_order, double span, size_t ndims) {
    this->ndims = ndims;
    loess->clear();
    label_order->clear();
    int find_index=0;
    int last_label=-10;
    int last_index=-10;
    int curr_label=0;
    //std::vector<PointCloud> labelled_points;  // 
    std::map<size_t, PointCloud> clusters;
    PointCloud c;
    for(size_t i = 0; i<labels->size(); i++) {
        curr_label=labels->at(i);
        if(curr_label!=-1) {
            if (clusters.find(curr_label) == clusters.end()) {
                clusters[curr_label] = PointCloud();
                label_order->push_back(curr_label);
            }
            clusters[curr_label].push_back(pc->at(i));
        }
    }

    /* if(smooth_radius>0) {
       for(auto &p : clusters) {
            c.clear();
            loess->push_back(c);
            smoothen_cluster(p.second, smooth_radius);
            do_loess(&p.second, &(loess->at(loess->size()-1)), 0.5);
        }
    } else { */
        for(auto &p : clusters) {
            if(p.first!=-1) {
                c.clear();
                loess->push_back(c);
                do_loess(&p.second, &(loess->at(loess->size()-1)), span);
            }
        }
    //}

}

void LOESS::do_loess(PointCloud * pc, PointCloud * out, double span) {
    size_t n = pc->size();
    double * x = new double[n];
    double * y = new double[n];
    double * t = new double[n];

    std::sort(pc->begin(), pc->end(), compareByT);

    double start=pc->at(0).t;
    double end=pc->at(pc->size()-1).t;

    for(size_t i = 0; i<n; i++) {
        x[i] = pc->at(i).data.at(0);
        y[i] = pc->at(i).data.at(1);
        t[i] = pc->at(i).t;
    }


    WeightedLowess::Options opt;

    opt.span=span;
    opt.anchors = pc->size();
    opt.num_threads=8;
    //size_t n_points=pc->size();
    double intervals = opt.anchors/(end-start);
    std::vector<double> fittedx(n), residsx(n);
    std::vector<double> fittedy(n), residsy(n);
    auto xwindows = WeightedLowess::define_windows(n-1, t, opt);
    WeightedLowess::compute(n, t, xwindows, y, fittedy.data(), residsy.data(), opt);
    WeightedLowess::compute(n, t, xwindows, x, fittedx.data(), residsx.data(), opt);

    if(ndims==4) {
        double * z = new double[n];
        for(size_t i = 0; i<n; i++) {
            z[i] = pc->at(i).data.at(2);
        }
        std::vector<double> fittedz(n), residsz(n);
        WeightedLowess::compute(n, t, xwindows, z, fittedz.data(), residsz.data(), opt);


        Point0 a=Point0(ndims);
        for(size_t i = 0; i<n; i++) {
            a.data.at(0) = fittedx.at(i);
            a.data.at(1) = fittedy.at(i);
            a.data.at(2) = fittedz.at(i);
            //a.t = pc->at(i).t;
            //a.t = xwindows.anchors.at(i);
            a.t = start+xwindows.anchors.at(i)*intervals;
            //a.t = start+i*intervals;
            out->push_back(a);
        }


        delete [] x;
        delete [] y;
        delete [] z;
        delete [] t;

    } else {
        Point0 a=Point0(3);
        for(size_t i = 0; i<(n-1); i++) {
            a.data.at(0) = fittedx.at(i);
            a.data.at(1) = fittedy.at(i);
            //a.t = pc->at(i).t;
            //a.t = xwindows.anchors.at(i);
            a.t = start+xwindows.anchors.at(i)*intervals;
            //a.t = start+i*intervals;
            out->push_back(a);
        }


        delete [] x;
        delete [] y;
        delete [] t;

    }

}