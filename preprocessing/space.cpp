#include "space.hpp"
#include <flann/flann.hpp>


Space::Space(PointCloud * pc, PointCloud * output, double eps) {
    this->eps = eps; 
    this->pc = pc;
    this->output = output;
    this->output->clear();
}


void Space::preprocess() {
    int n = (int)pc->size();
    if(!pseudorandom_order) {
        bool * to_be_kept = new bool[n];
        to_be_kept[0] = true;
        for(int i = 1; i<n; i++) {
            to_be_kept[i] = false;
        }
        // build indexed KDTREE here
        // first write data into helper array
        float * data = new float[n*3];
        int temp;
        for(int i=0; i<n; i++) {
            temp=i*3;
            *((data + temp)) = (float)pc->at(i).x;
            *((data + temp)+1) = (float)pc->at(i).y;
            *((data + temp)+2) = (float)pc->at(i).t;
        }

        flann::Matrix<float> dataset=flann::Matrix<float>(data, n, 3);   // construct dataset 
        flann::Matrix<float> query=flann::Matrix<float>(data, n, 3);     // construct queries (all datapoints again)

        flann::Index<flann::L2<float> > index(dataset, flann::KDTreeIndexParams(4));  
        index.buildIndex();                                             // build KdTree with Index 
        // initialize datastructures for radiusSearch
        std::vector<std::vector<size_t>> indices;
        std::vector<std::vector<float>> dists;

        // radius search for all points in pc
        index.radiusSearch(query, indices, dists, (float)eps, flann::SearchParams(128));

        output->push_back(pc->at(0));
        bool isolated=true;
        for(int i = 1; i<n; i++)  {
            isolated=true;
            for(int j = 0; j<(int)(indices.at(i).size()); j++) {
                if (to_be_kept[indices.at(i).at(j)]) {
                    isolated=false;
                    break;
                }
            }
            to_be_kept[i]=isolated;
            if(isolated) output->push_back(pc->at(i));
        }

    } 
} 