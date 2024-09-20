#include "assign.hpp"

#include <flann/flann.hpp>

Assign::Assign(PointCloud * downsampled, PointCloud * full, std::vector<int> * labels_downsampled, std::vector<int> *labels_full, size_t ndims) {

    size_t n_downsampled = downsampled->size();
    float * data_downsampled = new float[n_downsampled*ndims];
    int temp;
    for(int i=0; i<n_downsampled; i++) {
        temp=i*ndims;
        for(size_t j = 0; j<(ndims-1); j++) {
            *((data_downsampled + temp)+j) = (float)downsampled->at(i).data.at(j);  
        }
        *((data_downsampled + temp)+(ndims-1)) = (float)downsampled->at(i).t;
    }

    flann::Matrix<float> dataset=flann::Matrix<float>(data_downsampled, n_downsampled, ndims);   // construct dataset 

    //flann::AutotunedIndexParams ip = flann::AutotunedIndexParams(0.9, 0,0,0.1); 

    flann::Index<flann::L2<float> > index(dataset, flann::KDTreeIndexParams(2));  
    //flann::Index<flann::L2<float> > index(dataset, flann::LinearIndexParams());  
    index.buildIndex();                                              // build KdTree with Index 


    size_t n_full = full->size();
    float * data_full = new float[n_full*ndims];
    for(int i = 0; i<n_full; i++) {
        temp=i*ndims;
        for(size_t j = 0; j<(ndims-1); j++) {
            *((data_full + temp)+j) = (float)full->at(i).data.at(j);  
        }
        *((data_full + temp)+(ndims-1)) = (float)full->at(i).t;
    }

    flann::SearchParams sp = flann::SearchParams(124);
    //flann::SearchParams sp = flann::SearchParams(32);
    sp.cores=0;


    flann::Matrix<float> query=flann::Matrix<float>(data_full, n_full, ndims);     // construct queries (all datapoints again)

    std::vector<std::vector<size_t>> indices;
    std::vector<std::vector<float>> dists;
    index.knnSearch(query, indices, dists, 1, sp);
    for(size_t i = 0; i<full->size(); i++) {
        labels_full->push_back(labels_downsampled->at(indices.at(i).at(0)));
    }

    delete [] data_downsampled;
    delete [] data_full;

}


Assign::Assign(PointCloud * downsampled, PointCloud * full, std::vector<int> * labels_downsampled, std::vector<int> *labels_full, std::vector<size_t>* noise_points, size_t ndims) {

    size_t n_downsampled = downsampled->size();
    float * data_downsampled = new float[n_downsampled*ndims];
    int temp;
    for(int i=0; i<n_downsampled; i++) {
        temp=i*ndims;
        for(size_t j = 0; j<(ndims-1); j++) {
            *((data_downsampled + temp)+j) = (float)downsampled->at(i).data.at(j);  
        }
        *((data_downsampled + temp)+(ndims-1)) = (float)downsampled->at(i).t;
    }

    flann::Matrix<float> dataset=flann::Matrix<float>(data_downsampled, n_downsampled, ndims);

    //flann::AutotunedIndexParams ip = flann::AutotunedIndexParams(0.9, 0,0,0.1); 

    flann::Index<flann::L2<float> > index(dataset, flann::KDTreeIndexParams(2));  
    //flann::Index<flann::L2<float> > index(dataset, flann::LinearIndexParams());  
    index.buildIndex();                                              // build KdTree with Index 


    size_t n_full = full->size();
    //float * data_full = new float[(n_full-noise_points->size())*3];
    float * data_full = new float[n_full*ndims];
    size_t curr_idx_noise=0;
    for(size_t i = 0; i<n_full; i++) {
        temp=i*ndims;
        for(size_t j = 0; j<(ndims-1); j++) {
            *((data_full + temp)+j) = (float)full->at(i).data.at(j);  
        }
        *((data_full + temp)+(ndims-1)) = (float)full->at(i).t;
        //if((curr_idx_noise < noise_points->size()) && (noise_points->at(curr_idx_noise) == i)) {
        //    curr_idx_noise++;
        //} else {
            //temp=i*ndims;
            //*((data_full + temp)) = (float)full->at(i).x;
            //*((data_full + temp)+1) = (float)full->at(i).y;
            //*((data_full + temp)+2) = (float)full->at(i).t;
            //*((data_full + temp)) = (float)full->at((i-curr_idx_noise)).x;
            //*((data_full + temp)+1) = (float)full->at((i-curr_idx_noise)).y;
            //*((data_full + temp)+2) = (float)full->at((i-curr_idx_noise)).t;
        //}
    }

    flann::SearchParams sp = flann::SearchParams(124);
    //flann::SearchParams sp = flann::SearchParams(32);
    sp.cores=0;


    flann::Matrix<float> query=flann::Matrix<float>(data_full, n_full, ndims);     // construct queries (all datapoints again)

    std::vector<std::vector<size_t>> indices;
    std::vector<std::vector<float>> dists;
    index.knnSearch(query, indices, dists, 1, sp);
    
    curr_idx_noise=0;
    for(size_t i = 0; i<full->size(); i++) {
        if((curr_idx_noise < noise_points->size()) && (noise_points->at(curr_idx_noise) == i)) {
            curr_idx_noise++;
            labels_full->push_back(-1);
        } else {
            labels_full->push_back(labels_downsampled->at(indices.at(i).at(0)));
        }
    }

    delete [] data_downsampled;
    delete [] data_full;

}