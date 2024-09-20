    #include "mst.hpp"
    #include <flann/flann.hpp>
    #include <sstream>
    #include <fstream>
    #include <limits>
    #include <deque>
    #include <unordered_map>
    #include <utility>
    #include "../utils/mathstuff.hpp"


    PointCloud n1;
    PointCloud n2;
    // für test_kreuzung2 mit octree target_size 10
    //int idx123=10;     // an dem Knick
    //int idx123=85;     // an einer anderen inkonsistenten Kante
    //int idx123=76;
    //int idx123=30;       // innerhalb einer langen Strecke 
    //int idx123=200;      // innerhalb einer kurzen Strecke
    //int idx123=10;
    //int idx123=144;        // innerhalb einer gebogenen Strecke
    //int idx123=101;
    //int idx123=6;            // super große Kosinus-Ähnlichkeit
    //int idx123=143;          // super kleine Kosinus-Ähnlichkeit
    int idx123=103;
    //int idx123=55;
    PointCloud pca_idx123_n1;
    PointCloud pca_idx123_n2; 


    void RANSAC(std::vector<Point> * neighbourhood, float t, size_t * i, size_t * j) {
        MathStuff ms;
        size_t max_i=0;
        size_t max_j=0;
        size_t max_num_supporters=0;
        size_t curr_supporters=0;
        Point a;
        Point b;
        Point c;
        for(size_t i = 0; i<neighbourhood->size(); i++) {
            a=neighbourhood->at(i);
            for(size_t j = 0; j<i; j++) {
                b=neighbourhood->at(j)-neighbourhood->at(i);
                curr_supporters=0;
                for(size_t k = 0; k<neighbourhood->size(); k++) {
                    c=neighbourhood->at(k);
                    //std::cout<<ms.orthogonal_lsq_distance(&a,&b,&c)<<"\n";
                    if(ms.orthogonal_lsq_distance(&a,&b,&c)<t) curr_supporters++;
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
                    if(ms.orthogonal_lsq_distance(&a,&b,&c)<t) curr_supporters++;
                }
                if(curr_supporters>max_num_supporters) {
                    max_i=i;
                    max_j=j;
                    max_num_supporters=curr_supporters;
                }
            }
        }
        *i=max_i;
        *j=max_j;
    }




    DirectedMST::DirectedMST(PointCloud * pc, std::vector<int> * labels, double max_edge_dist, size_t max_branch_depth, bool edge_length_criterion, bool pca_angle_criterion, double pca_angle_value, int octree_target_size) {

        labels->clear();
        for(size_t i = 0; i<pc->size(); i++) {
            labels->push_back(-1);
            mst.push_back(DirectedNode(i));
        }
        this->edge_length_criterion=edge_length_criterion;
        this->pca_angle_criterion=pca_angle_criterion;
        this->twoD = twoD;
        this->max_edge_dist=max_edge_dist;
        this->max_branch_depth=max_branch_depth;
        this->pca_angle_value=pca_angle_value;
        this->labels = labels;
        this->pc = pc;
        this->octree_target_size=octree_target_size;
        build_MST();
        //remove_time_gaps();
        remove_long_edges();
        remove_knees();

        //remove_branches();

        color_nodes();
        const char * fname = "pca_plot.gnuplot";
        draw_pca(fname);
        const char * fname2 = "plot_edges.gnuplot";
        draw_pca_debug(fname2, false);
        //const char * fname3 = "pca_plot_edges2D.gnuplot";
        //draw_pca_debug(fname3, true);
    }


    void DirectedMST::build_MST() { // NO FLANN
        size_t n = pc->size();
        Point curr;
        const size_t big_size_t=std::numeric_limits<size_t>::max();
        const double big_double=std::numeric_limits<double>::max();
        double min_distance=big_double;
        size_t min_index=big_size_t;
        double curr_distance=0;
        double delta_t_squared=0;

        for(size_t i = 0; i<n; i++) {
            curr = pc->at(i);
            min_distance=big_double;
            min_index=big_size_t;
            for(size_t j = i+1; j<n; j++) {
                delta_t_squared = (curr.t-pc->at(j).t)*(curr.t-pc->at(j).t);
                if(delta_t_squared > min_distance) break;
                curr_distance=curr.no_root_euclidian_distance(pc->at(j));
                if(min_distance > curr_distance) {
                    min_distance=curr_distance;
                    min_index=j;
                }
            }
            if(min_distance!=big_double) {
                mst.at(min_index).predecessor.push_back(&(mst.at(i)));
                mst.at(min_index).dist_predecessor.push_back(min_distance);
                mst.at(i).successor = &(mst.at(min_index));
                mst.at(i).dist_successor = min_distance;
            }
        }

    }


    void DirectedMST::remove_time_gaps() {
        for(size_t i = 0; i<mst.size(); i++) {
            if(mst.at(i).successor) {
                if(std::abs(pc->at(i).t - pc->at(mst.at(i).successor->index).t)>50) {
                    delete_edge(i);
                }
            }
        }
    }


    void DirectedMST::remove_long_edges() {
        std::vector<size_t> snipsnap;  
        size_t a, b;
        double global_max_dist_threshold = choose_global_max_dist_threshold();
        double local_max_dist_threshold=global_max_dist_threshold;
        std::cout<<global_max_dist_threshold<<"\n";
        //std::cout<<"Chose global max_dist_threshold="<<std::sqrt(global_max_dist_threshold)<<"\n";
        for(size_t i = 0; i<mst.size(); i++) {
            //if(this->max_edge_dist==-2) local_max_dist_threshold=choose_local_max_dist_threshold(i);
            if(mst.at(i).successor) {
                if((mst.at(i).dist_successor > local_max_dist_threshold)) {
                    //std::cout<<"SNIP SNAP\n";
                    //snipsnap.push_back(i);
                    delete_edge(i);
                }
            }
        }

        /* for(size_t i = 0; i<mst.size(); i++) {
            if(mst.at(i).successor) {
                if(this->max_edge_dist==-2) local_max_dist_threshold=choose_local_max_dist_threshold(i);
                if((mst.at(i).dist_successor > local_max_dist_threshold)) {
                    delete_edge(i);
                }
            }
        } */

        /* for(size_t i : snipsnap) {
            delete_edge(i);
        } */
    }



    void DirectedMST::remove_knees() {
        const float ransac_t=300.0;
        //std::vector<size_t> to_be_removed;
        const size_t max_depth=5;
        //size_t curr_depth=0;
        PointCloud vectors_n1;
        PointCloud vectors_n2;
        //std::deque<DirectedNode*> open;
        //std::deque<size_t> depth;
        //std::unordered_map<size_t, bool> visited;
        std::vector<Point> ps;
        //DirectedNode * curr_node;
        Point a_n1;
        Point b_n1;
        Point a_n2;
        Point b_n2;
        MathStuff mathstuff;
        double cos_similarity_n1=0.0;
        double cos_similarity_n2=0.0;
        for(size_t index = 0; index<mst.size(); index++) {
            if(mst.at(index).successor && !((pc->at(index)-pc->at(mst.at(index).successor->index)) == Point(0,0,0))) {
                /* curr_depth=0;
                open.clear();
                depth.clear();
                open.push_back(mst.at(index).successor);
                depth.push_back(0);

                // first the next point
                vectors_n1.clear();
                while(open.size()>0) {
                    curr_node = open.at(0);
                    open.pop_front();
                    curr_depth = depth.at(0);
                    depth.pop_front();
                    if( (curr_depth < max_depth) && curr_node->successor ) {
                        vectors_n1.push_back(pc->at(curr_node->index));
                        open.push_back(curr_node->successor);
                        depth.push_back(curr_depth+1);
                    }
                }
                // then all previous points
                vectors_n2.clear();
                for(size_t i = 0; i<mst.at(index).predecessor.size(); i++) {
                    open.push_back(mst.at(index).predecessor.at(i));
                    depth.push_back(0);
                }
                while(open.size()>0) {
                    curr_node = open.at(0);
                    open.pop_front();
                    curr_depth = depth.at(0);
                    depth.pop_front();
                    if( (curr_depth < max_depth)  ) {
                        for(size_t i = 0; i<curr_node->predecessor.size(); i++) {
                            vectors_n2.push_back(pc->at(curr_node->index));
                            open.push_back(curr_node->predecessor.at(i));
                            depth.push_back(curr_depth+1);
                        }
                    }
                } */

                BFS(index, &vectors_n1, &vectors_n2, max_depth);


                /* int num_doubles_n1=0;
                for(size_t i = 0; i<vectors_n1.size(); i++) {
                    for(size_t j = 0; j<i; j++) {
                        if(vectors_n1.at(i)==vectors_n1.at(j)) num_doubles_n1++;
                    }
                    for(size_t j = i+1; j<vectors_n1.size(); j++) {
                        if(vectors_n1.at(i)==vectors_n1.at(j)) num_doubles_n1++;
                    }
                }

                int num_doubles_n2=0;
                for(size_t i = 0; i<vectors_n2.size(); i++) {
                    for(size_t j = 0; j<i; j++) {
                        if(vectors_n2.at(i)==vectors_n2.at(j)) num_doubles_n2++;
                    }
                    for(size_t j = i+1; j<vectors_n2.size(); j++) {
                        if(vectors_n2.at(i)==vectors_n2.at(j)) num_doubles_n2++;
                    }
                }
                std::cout<<num_doubles_n1<<", "<<num_doubles_n2<<"\n"; */



                if((vectors_n1.size()>=(size_t)(max_depth/2)) && (vectors_n2.size()>=(size_t)(max_depth/2))) {
                    size_t idx_i;
                    size_t idx_j;
                    RANSAC(&vectors_n1, ransac_t, &idx_i, &idx_j);
                    a_n1 = pc->at(index);
                    b_n1 = vectors_n1.at(idx_i)-vectors_n1.at(idx_j);
                    RANSAC(&vectors_n2, ransac_t, &idx_i, &idx_j);
                    a_n2 = pc->at(mst.at(index).successor->index);
                    b_n2 = vectors_n2.at(idx_i)-vectors_n2.at(idx_j);
                    /* mathstuff.orthogonal_lsq(vectors_n1, &a_n1, &b_n1);
                    mathstuff.orthogonal_lsq(vectors_n2, &a_n2, &b_n2); */
                    ps.clear();
                    ps.push_back(a_n1);
                    ps.push_back(b_n1);
                    pca_equations_n1.push_back(ps);
                    ps.clear();
                    ps.push_back(a_n2);
                    ps.push_back(b_n2);
                    pca_equations_n2.push_back(ps);
                    cos_similarity_n1 = (pc->at(index)-pc->at(mst.at(index).successor->index)).cosine_similarity(b_n1);
                    cos_similarity_n2 = (pc->at(index)-pc->at(mst.at(index).successor->index)).cosine_similarity(b_n2);
                    kneeness.push_back(std::min(std::abs(cos_similarity_n1), std::abs(cos_similarity_n2))); 
                    if(std::min(std::abs(cos_similarity_n1), std::abs(cos_similarity_n2)) < pca_angle_value) {
                        delete_edge(index);
                        //std::cout<<index<<"\n";
                        //std::cout<<cos_similarity_n1<<" : "<<cos_similarity_n2<<"\n";
                    } /* else {
                        std::cout<<"hier\n";
                        std::cout<<index<<"\n";
                        std::cout<<cos_similarity_n1<<" : "<<cos_similarity_n2<<"\n";
                    } */
                } else if(vectors_n1.size()>=(size_t)(max_depth/2)) {
                    //mathstuff.orthogonal_lsq(vectors_n1, &a_n1, &b_n1);
                    size_t idx_i;
                    size_t idx_j;
                    RANSAC(&vectors_n1, ransac_t, &idx_i, &idx_j);
                    a_n1 = vectors_n1.at(idx_i);
                    b_n1 = vectors_n1.at(idx_i)-vectors_n1.at(idx_j);
                    ps.clear();
                    ps.push_back(a_n1);
                    ps.push_back(b_n1);
                    pca_equations_n1.push_back(ps);
                    cos_similarity_n1 = (pc->at(index)-pc->at(mst.at(index).successor->index)).cosine_similarity(b_n1);
                    kneeness.push_back(std::abs(cos_similarity_n1)); 
                    if(std::abs(cos_similarity_n1) < pca_angle_value) {
                        delete_edge(index);
                        //std::cout<<index<<"\n";
                        //std::cout<<cos_similarity_n1<<"\n";
                    } /* else {
                        std::cout<<"hier\n";
                        std::cout<<index<<"\n";
                        std::cout<<cos_similarity_n1<<"\n";
                    } */
                } else if(vectors_n2.size()>=(size_t)(max_depth/2)) {
                    //mathstuff.orthogonal_lsq(vectors_n2, &a_n2, &b_n2);
                    size_t idx_i;
                    size_t idx_j;
                    RANSAC(&vectors_n2, ransac_t, &idx_i, &idx_j);
                    a_n2 = vectors_n2.at(idx_i);
                    b_n2 = vectors_n2.at(idx_i)-vectors_n2.at(idx_j);
                    ps.clear();
                    ps.push_back(a_n2);
                    ps.push_back(b_n2);
                    pca_equations_n1.push_back(ps);
                    cos_similarity_n2 = (pc->at(index)-pc->at(mst.at(index).successor->index)).cosine_similarity(b_n2);
                    kneeness.push_back(std::abs(cos_similarity_n2)); 
                    if(std::abs(cos_similarity_n2) < pca_angle_value) {
                        delete_edge(index);
                        //std::cout<<index<<"\n";
                        //std::cout<<cos_similarity_n2<<"\n";
                    } /* else {
                        std::cout<<"hier\n";
                        std::cout<<index<<"\n";
                        std::cout<<cos_similarity_n2<<"\n";
                    } */
                } else {
                    kneeness.push_back(-2.0);
                }

                if(index==idx123) {
                    n1 = vectors_n1;
                    n2 = vectors_n2;
                    if(kneeness.at(kneeness.size()-1) > 0) {
                        pca_idx123_n1.push_back(a_n1);
                        pca_idx123_n1.push_back(b_n1);
                        pca_idx123_n2.push_back(a_n2);
                        pca_idx123_n2.push_back(b_n2);
                    }
                }

            } else {
                kneeness.push_back(-1.0);
            }
        }
        //for(size_t index : to_be_removed) {
        //    delete_edge(index);
        //}
    }       


    ssize_t DirectedMST::BFS(size_t index, PointCloud * n1, PointCloud * n2, size_t max_depth) {
        n1->clear();
        n2->clear();
        size_t curr_depth=0;
        DirectedNode * curr_node;
        std::deque<DirectedNode*> open;
        std::deque<size_t> depth;
        std::unordered_map<size_t, bool> closed;
        // !!!!!!!!!!!!!!!!!! n1 !!!!!!!!!!!!!!!!!!!
        for(size_t i = 0; i<mst.at(index).predecessor.size(); i++) {
            //n1->push_back(pc->at(index));
            open.push_back(&mst.at(index));
            depth.push_back(0);
        }

        while(!open.empty()) {
            curr_node = open.at(0);
            open.pop_front();
            curr_depth = depth.at(0);
            depth.pop_front();
            if( (curr_depth < max_depth) ) {
                n1->push_back(pc->at(curr_node->index));
                for(size_t i = 0; i<curr_node->predecessor.size(); i++) {
                    /* for(size_t j = 0; j<n1->size(); j++) {
                        if(n1->at(j)==pc->at(curr_node->predecessor.at(i)->index)) {
                            std::cout<<"--------------\n";
                            std::cout<<j<<"\n";
                            std::cout<<n1->size()<<"\n";
                            std::cout<<curr_node->predecessor.at(i)->index<<"\n";
                            std::cout<<curr_node->index<<"\n";
                            break;
                        }
                    } */
                    
                    //for(size_t j = 0; j<curr_node->predecessor.at(i)->predecessor.size(); j++) {
                    open.push_back(curr_node->predecessor.at(i));
                    depth.push_back(curr_depth+1);
                    //}
                }
            }
        }

        /* if(index==11) {
            for(size_t i = 0; i<n1->size(); i++) {
                std::cout<<n1->at(i)<<"\n";
            }
        } */


        /* int num_doubles_n1=0;
        for(size_t i = 0; i<n1->size(); i++) {
            for(size_t j = 0; j<i; j++) {
                if(n1->at(i)==n1->at(j)) num_doubles_n1++;
            }
            for(size_t j = i+1; j<n1->size(); j++) {
                if(n1->at(i)==n1->at(j)) num_doubles_n1++;
            }
        }
        std::cout<<num_doubles_n1<<"\n"; */


        // !!!!!!!!!!!!!!!!!! n2 !!!!!!!!!!!!!!!!!!!

        closed[mst.at(index).index] = true;
        open.push_back(mst.at(index).successor);
        depth.push_back(0);
        n2->push_back(pc->at(mst.at(index).successor->index));
        while(!open.empty()) {
            curr_node = open.at(0);
            open.pop_front();
            curr_depth = depth.at(0);
            depth.pop_front();

            if( curr_depth < max_depth ) {
                if(curr_node->successor && (closed.find(curr_node->successor->index) == closed.end())) {
                    n2->push_back(pc->at(curr_node->successor->index));
                    open.push_back(curr_node->successor);
                    depth.push_back(curr_depth+1);
                }
                /* for(size_t i = 0; i<curr_node->predecessor.size(); i++) {
                    if(closed.find(curr_node->predecessor.at(i)->index) == closed.end()) {
                        n2->push_back(pc->at(curr_node->predecessor.at(i)->index));
                        open.push_back(curr_node->predecessor.at(i));
                        depth.push_back(curr_depth+1);
                    }
                } */
            }
            closed[curr_node->index] = true;
        }

        /* if(index==11) {
            std::cout<<"---------------\n";
            for(size_t i = 0; i<n2->size(); i++) {
                std::cout<<n2->at(i)<<"\n";
            }
        } */


        return n1->size()+n2->size();
    }


    void DirectedMST::delete_edge(size_t index) {
        DirectedNode * a = &mst.at(index);

        removed_edges.push_back(pc->at(a->index));
        if(mst.at(index).successor) {
            DirectedNode * b = mst.at(index).successor;
            removed_edges.push_back(pc->at(b->index));
            std::vector<DirectedNode *> temp_predecessor=b->predecessor;
            b->predecessor.clear();
            for(size_t i = 0; i<temp_predecessor.size(); i++) {
                if(temp_predecessor.at(i)!=a) b->predecessor.push_back(temp_predecessor.at(i));
            }
        }
        a->successor = NULL;
    }

    double DirectedMST::choose_local_max_dist_threshold(size_t index) { // index is the index of the node where the edge in question is starting
        const size_t max_depth=10;  // die Hardgecodedness muss noch verändert werden
        const double c=10;
        size_t curr_depth=0;
        std::vector<double> distances_n1;
        std::vector<double> distances_n2;
        std::deque<DirectedNode*> open;
        std::deque<size_t> depth;
        DirectedNode * curr_node;
        if(mst.at(index).successor) { 
            open.push_back(mst.at(index).successor);
            depth.push_back(0);
        }

        // first the next point
        while(open.size()>0) {
            curr_node = open.at(0);
            open.pop_front();
            curr_depth = depth.at(0);
            depth.pop_front();
            if( (curr_depth < max_depth) && curr_node->successor ) {
                distances_n1.push_back(curr_node->dist_successor);
                open.push_back(curr_node->successor);
                depth.push_back(curr_depth+1);
            }
        }
        // then all previous points
        for(size_t i = 0; i<mst.at(index).predecessor.size(); i++) {
            open.push_back(mst.at(index).predecessor.at(i));
            depth.push_back(0);
        }
        while(open.size()>0) {
            curr_node = open.at(0);
            open.pop_front();
            curr_depth = depth.at(0);
            depth.pop_front();
            if( (curr_depth < max_depth)  ) {
                for(size_t i = 0; i<curr_node->predecessor.size(); i++) {
                    distances_n2.push_back(curr_node->dist_predecessor.at(i));
                    open.push_back(curr_node->predecessor.at(i));
                    depth.push_back(curr_depth+1);
                }
            }
        }

        double mean1=0.0;
        for(size_t i = 0; i<distances_n1.size(); i++) {
            mean1+=distances_n1.at(i);
        }
        if(distances_n1.size() > 0) mean1/=distances_n1.size();
        else mean1=0;
        double sd1=0.0;
        for(size_t i = 0; i<distances_n1.size(); i++) {
            sd1+=(distances_n1.at(i)-mean1)*(distances_n1.at(i)-mean1);
        }
        if(distances_n1.size() > 1) sd1=std::sqrt(sd1/(distances_n1.size()-1));
        else sd1=0;


        double mean2=0.0;
        for(size_t i = 0; i<distances_n2.size(); i++) {
            mean2+=distances_n2.at(i);
        }
        if(distances_n2.size()>0) mean2/=distances_n2.size();
        else mean2 = 0;
        double sd2=0.0;
        for(size_t i = 0; i<distances_n2.size(); i++) {
            sd2+=(distances_n2.at(i)-mean2)*(distances_n2.at(i)-mean2);
        }
        if(distances_n2.size()>1) sd2=std::sqrt(sd2/(distances_n2.size()-1));
        else sd2=0;

        //if(distances_n1.size() > distances_n2.size()) return mean1+sd1;
        //else return mean2+sd2;
        //std::cout<<"----------\n";
        //std::cout<<distances_n1.size()<<", "<<distances_n2.size()<<"\n";
        //std::cout<<mean1+sd1<<", "<<mean2+sd2<<"\n";

        //if(mean1==0) mean1 = std::numeric_limits<double>::max();
        //if(mean2==0) mean2 = std::numeric_limits<double>::max();
        //return std::min(mean1+c*sd1, mean2+c*sd2);

        std::sort( distances_n1.begin(), distances_n1.end());
        std::sort( distances_n2.begin(), distances_n2.end());
        double median1=0.0;
        double median2=0.0;
        if(distances_n1.size()>0) median1 = distances_n1.at((int)(distances_n1.size()/2));
        if(distances_n2.size()>0) median2 = distances_n2.at((int)(distances_n2.size()/2));
        return std::max(median1+c*sd1, median2+c*sd2);

        return std::max(mean1+c*sd1, mean2+c*sd2);
    }


    double DirectedMST::choose_global_max_dist_threshold() {
        if(octree_target_size>0) {
            double val=4*std::sqrt(2*(double)(octree_target_size*octree_target_size));
            return val*val;
        }
        const double c=1;
        if(this->max_edge_dist >= 0)  return this->max_edge_dist;
        else {
            double mean=0.0;
            size_t count=0.0;
            for(size_t i = 0; i<mst.size(); i++) {
                if(mst.at(i).successor) {
                    mean+=mst.at(i).dist_successor;
                    count++;
                }
            }
            if(count>0) mean/=count;
            else std::cout<<"Error: choose_max_dist_threshold called before a tree was build\n";
            double sd=0.0;
            for(size_t i = 0; i<mst.size(); i++) {
                if(mst.at(i).successor) sd+=(mst.at(i).dist_successor-mean)*(mst.at(i).dist_successor-mean);
            }
            sd=std::sqrt(sd/(count-1));
            return mean+c*sd;
        }
    }


    /* double MST::calc_max_edge_dist_estimate(std::vector<std::vector<float>> * dists) {
        std::vector<float> distances;
        for(size_t i = 0; i<dists->size(); i++) {
            for(size_t j = 1; j<dists->at(i).size(); j++) {
                distances.push_back(dists->at(i).at(j));
            }
        }
        std::sort( distances.begin(), distances.end());
        float last_difference=distances.at(1)-distances.at(0);
        float difference=0.0;
        
        //for(size_t i = 2; i<distances.size(); i++) {
        //    difference=distances.at(i) - distances.at(i-1);
        //    if((difference > (10000*last_difference)) && (last_difference>0)) {
        //        return distances.at(i);
        //    }
        //    last_difference=difference;
        //}
        return (double)distances.at(std::round(distances.size()*0.95));
    } */


    void DirectedMST::calc_pca(std::vector<std::vector<size_t>> * indices) {
        MathStuff mathstuff;
        Point a;
        Point b;
        PointCloud cloud;
        Point curr_point;
        std::vector<double> angles;
        std::vector<Point> ps;
        for(size_t i = 0; i<indices->size(); i++) {
            curr_point=pc->at(i);
            cloud.clear();
            for(size_t j = 0; j<indices->at(i).size(); j++) {
                cloud.push_back(pc->at(indices->at(i).at(j)));
            }
            mathstuff.orthogonal_lsq(cloud, &a, &b);
            ps.clear();
            ps.push_back(a);
            ps.push_back(b);
            pca_equations.push_back(ps);
        }
    }


    double DirectedMST::inconsistent_edge_length(size_t index_a , size_t index_b, std::vector<float> * dists_a, std::vector<float> * dists_b, std::vector<size_t> * indices_a, std::vector<size_t> * indices_b) {
        std::vector<float> dists;
        size_t i=0;
        size_t j=0;
        double mean=0.0;
        while((i<dists_a->size()) && (j<dists_b->size())) {
            if(dists_a->at(i) <= dists_b->at(j)) {
                //if(indices_a->at(i) > index_a) {
                    mean+=dists_a->at(i);
                    dists.push_back(dists_a->at(i));
                //}
                i++;
            } else {
                //if(indices_b->at(j) > index_b) {
                    mean+=dists_b->at(j);
                    dists.push_back(dists_b->at(j));
                //}
                j++;
            }
        }
        while(i<dists_a->size()) {
            mean+=dists_a->at(i);
            dists.push_back(dists_a->at(i));
            i++;
        }
        while(j<dists_b->size()) {
            mean+=dists_b->at(j);
            dists.push_back(dists_b->at(j));
            j++;
        }
        //std::cout<<"dists.size="<<dists.size()<<"\n";

        double median=dists.at((int)(dists.size()*0.5));

        mean=mean/(dists.size());
        double iqr=dists.at((int)(dists.size()*0.75)) - dists.at((int)(dists.size()*0.25));
        return median+iqr;
    }


    void DirectedMST::remove_branches() {
        bool remove_edge=false;
        size_t n = mst.size();
        std::vector<size_t> to_be_removed; 
        double min_distance=0.0;
        DirectedNode * min_ptr=NULL;
        int removed_edges=0;
        size_t big_number=std::numeric_limits<size_t>::max();
        for(size_t i = 0; i<n; i++) {
            to_be_removed.clear();
            if(mst.at(i).predecessor.size() >= 2) {
                for(size_t j = 0; j < mst.at(i).predecessor.size(); j++) {
                    remove_edge=DFS(mst.at(i).predecessor.at(j),0);
                    if(remove_edge) {
                        to_be_removed.push_back(j);
                        removed_edges++;
                    }
                }
                min_distance=big_number;
                for(size_t j : to_be_removed) {
                    //if(mst.at(i).dist_predecessor.at(j) < min_distance) {
                    //    min_distance = mst.at(i).dist_predecessor.at(j);
                    //    min_ptr = mst.at(i).predecessor.at(j);
                    //}
                    mst.at(i).predecessor.at(j)->successor = NULL;
                }
                //if(min_ptr) {
                //    mst.at(i).predecessor.clear();
                //    mst.at(i).predecessor.push_back(min_ptr);
                //    min_ptr->successor = &mst.at(i);
                //}
            }
        } 
        //std::cout<<"Removed: "<<removed_edges<<" edges\n";
    }

    bool DirectedMST::DFS(DirectedNode * n, size_t curr_depth) {
        bool result=false;
        if(max_branch_depth >= curr_depth) {
            for(std::vector<DirectedNode *>::iterator it = n->predecessor.begin(); it!=n->predecessor.end(); ++it) {
                if(n->successor) 
                    result = DFS(n->successor, ++curr_depth);
                if(result) return result;
            }
            return result;
        } else {
            return true;
        }
    }


    int DirectedMST::expand(DirectedNode * node, int cluster_id) {
        if(labels->at(node->index) != -1) {
            return labels->at(node->index);
        }
        int collision=cluster_id;
        if(node->successor)
            collision = expand(node->successor, cluster_id);   
        labels->at(node->index) = collision;
        return collision;
    }


    void DirectedMST::color_nodes() {
        size_t n = pc->size();
        int cluster_id = 0;
        int collision=0;
        for(size_t i = 0; i<n; i++) {
            if(labels->at(i)==-1) {
                collision = expand(&mst.at(i), cluster_id);
                if(cluster_id==collision) cluster_id++; 
            } 
        }


        ////// labels with small amounts of data are noise?
    }


    /* bool MST::draw_mst(const char * fname) {
        // saves cloud in gnuplot file
        std::ofstream of(fname);
        of << std::fixed;  // set float style
        if (!of.is_open()) {
        std::cerr << "[Error] could not save under '" << fname << std::endl;
        return false;
        }

        of << "splot '-' using 1:2:3 with points lc 'black' pointtype 7 ps 1, '-' with lines palette title 'edges'; pause -1;\n";
        for(size_t i = 0; i<pc->size(); i++) {
            of<<pc->at(i).x<<" "<<pc->at(i).y<<" "<<pc->at(i).t<<"\n";
        }
        Point curr;
        Point successor;
        of <<"e\n";
        for(size_t i = 0; i<mst.size(); i++) { 
            curr = pc->at(i);
            if(mst.at(i).successor) {
                successor = pc->at(mst.at(i).successor->index);
                of<<curr.x<<" "<<curr.y<<" "<<curr.t<<" "<<mst.at(i).dist_successor<<"\n";
                of<<successor.x<<" "<<successor.y<<" "<<successor.t<<" "<<mst.at(i).dist_successor<<"\n\n\n";
            }
        }
        of.close();
        return true;
    }  */

    bool DirectedMST::draw_mst(const char * fname) {
        size_t version=0;
        // saves cloud in gnuplot file
        std::ofstream of(fname);
        of << std::fixed;  // set float style
        if (!of.is_open()) {
        std::cerr << "[Error] could not save under '" << fname << std::endl;
        return false;
        }

        if(version==0) of << "splot '-' using 1:2:3 with points lc 'black' pointtype 7 ps 1, '-' with lines palette title 'edges',  '-' with points lc 'blue', '-' with points lc 'green', '-' with points lc 'red', '-' with lines lc 'blue' title 'pca_n1', '-' with lines lc 'green' title 'pca_n2', '-' with lines lc 'red' title 'edge direction'; pause -1;\n";
        if(version==1) of << "set cbrange [0:1]\nsplot '-' using 1:2:3:4 with points palette pointtype 7 ps 1, '-' with lines title 'edges', '-' with points lc 'green' pointtype 5 pointsize 1; pause -1;\n";
        for(size_t i = 0; i<pc->size(); i++) {
            of<<pc->at(i).x<<" "<<pc->at(i).y<<" "<<pc->at(i).t<<" "<<kneeness.at(i)<<"\n";
        }
        Point curr;
        Point successor;
        of <<"e\n";
        for(size_t i = 0; i<mst.size(); i++) { 
            curr = pc->at(i);
            if(mst.at(i).successor) {
                successor = pc->at(mst.at(i).successor->index);
                of<<curr.x<<" "<<curr.y<<" "<<curr.t<<" "<<mst.at(i).dist_successor<<"\n";
                of<<successor.x<<" "<<successor.y<<" "<<successor.t<<" "<<mst.at(i).dist_successor<<"\n\n\n";
            }
        }
        of <<"e\n";
        
        if(version==1) {
            for(size_t i = 0; i<kneeness.size(); i++) {
                if(kneeness.at(i)<0) {
                    of<<pc->at(i).x<<" "<<pc->at(i).y<<" "<<pc->at(i).t<<"\n";
                    //std::cout<<i<<"\n";
                }
            }
        }

        if(version==0) { 
            for(size_t i = 0; i<n1.size(); i++) {
                of<<n1.at(i).x<<" "<<n1.at(i).y<<" "<<n1.at(i).t<<"\n";
            }
            of<<"e\n";
            for(size_t i = 0; i<n2.size(); i++) {
                of<<n2.at(i).x<<" "<<n2.at(i).y<<" "<<n2.at(i).t<<"\n";
            }
            of<<"e\n";
            of<<pc->at(idx123).x<<" "<<pc->at(idx123).y<<" "<<pc->at(idx123).t<<"\ne\n";

            float factor=20;

            if(pca_idx123_n1.size()==2) {
                of<<pca_idx123_n1.at(0).x<<" "<<pca_idx123_n1.at(0).y<<" "<<pca_idx123_n1.at(0).t<<"\n";
                of<<pca_idx123_n1.at(0).x+factor*pca_idx123_n1.at(1).x<<" "<<pca_idx123_n1.at(0).y+factor*pca_idx123_n1.at(1).y<<" "<<pca_idx123_n1.at(0).t+factor*pca_idx123_n1.at(1).t<<"\n\n\n";
            }
            of<<"e\n";
            if(pca_idx123_n2.size()==2) {
                of<<pca_idx123_n2.at(0).x<<" "<<pca_idx123_n2.at(0).y<<" "<<pca_idx123_n2.at(0).t<<"\n";
                of<<pca_idx123_n2.at(0).x+factor*pca_idx123_n2.at(1).x<<" "<<pca_idx123_n2.at(0).y+factor*pca_idx123_n2.at(1).y<<" "<<pca_idx123_n2.at(0).t+factor*pca_idx123_n2.at(1).t<<"\n\n\n";
            }
            of<<"e\n";
            if(mst.at(idx123).successor) {
                of<<pc->at(mst.at(idx123).index).x<<" "<<pc->at(mst.at(idx123).index).y<<" "<<pc->at(mst.at(idx123).index).t<<"\n";
                //of<<pc->at(mst.at(idx123).successor->index).x<<" "<<pc->at(mst.at(idx123).successor->index).y<<" "<<pc->at(mst.at(idx123).successor->index).t<<"\n";
                Point edge = pc->at(mst.at(idx123).index) - pc->at(mst.at(idx123).successor->index);
                of<<pc->at(mst.at(idx123).index).x + 20*edge.x<<" "<<pc->at(mst.at(idx123).index).y + 20*edge.y<<" "<<pc->at(mst.at(idx123).index).t + 20*edge.t<<"\n";
            }
            of<<"e\n";
        }
        of.close();
        return true;
    }    


    bool DirectedMST::draw_pca(const char * fname) {
        if(!pca_angle_criterion) {
            std::cout<<"Can't draw a pca plot without having calculated pca's\n";
            return false;
        }
        // saves cloud in gnuplot file
        std::ofstream of(fname);
        of << std::fixed;  // set float style
        if (!of.is_open()) {
            std::cerr << "[Error] could not save under '" << fname << std::endl;
            return false;
        }

        of << "splot '-' using 1:2:3 with points lc 'black' pointtype 7 ps 1, '-' with lines lc 'green' title 'n1', '-' with lines lc 'red' title 'n2'; pause -1;\n";
        for(size_t i = 0; i<pc->size(); i++) {
            of<<pc->at(i).x<<" "<<pc->at(i).y<<" "<<pc->at(i).t<<"\n";
        }
        Point a;
        Point b;
        of <<"e\n";
        double scale = 10.0;
        //std::cout<<pca_equations_n1.size()<<", "<<pca_equations_n2.size()<<"\n";
        for(size_t i = 0; i<pca_equations_n1.size(); i++) { 
            a = pca_equations_n1.at(i).at(0);
            b = pca_equations_n1.at(i).at(1);
            //std::cout<<a<<"\n";
            //of<<pc->at(i).x<<" "<<pc->at(i).y<<" "<<pc->at(i).t<<"\n";
            //of<<pc->at(i).x+scale*b.x<<" "<<pc->at(i).y+scale*b.y<<" "<<pc->at(i).t+scale*b.t<<"\n\n\n";
            of<<a.x<<" "<<a.y<<" "<<a.t<<"\n";
            of<<a.x+scale*b.x<<" "<<a.y+scale*b.y<<" "<<a.t+scale*b.t<<"\n\n\n";
        }
        
        of <<"e\n";
        scale = 10.0;
        for(size_t i = 0; i<pca_equations_n2.size(); i++) { 
            a = pca_equations_n2.at(i).at(0);
            b = pca_equations_n2.at(i).at(1);
            //std::cout<<a<<"\n";
            //of<<pc->at(i).x<<" "<<pc->at(i).y<<" "<<pc->at(i).t<<"\n";
            //of<<pc->at(i).x+scale*b.x<<" "<<pc->at(i).y+scale*b.y<<" "<<pc->at(i).t+scale*b.t<<"\n\n\n";
            of<<a.x<<" "<<a.y<<" "<<a.t<<"\n";
            of<<a.x+scale*b.x<<" "<<a.y+scale*b.y<<" "<<a.t+scale*b.t<<"\n\n\n";
        }
        of.close();
        return true;
    } 


    bool DirectedMST::draw_pca_debug(const char * fname, bool twoD) {
        // saves cloud in gnuplot file
        std::ofstream of(fname);
        of << std::fixed;  // set float style
        if (!of.is_open()) {
        std::cerr << "[Error] could not save under '" << fname << std::endl;
        return false;
        }

        if(twoD) {
            of << "plot '-' using 1:2 with points pointtype 7 ps 1, '-' with lines lc 'black' title 'edges', '-' with lines lc 'green' title 'removed_edges'; pause -1;\n";
            for(size_t i = 0; i<pc->size(); i++) {
                //of<<pc->at(i).x<<" "<<pc->at(i).y<<" "<<pc->at(i).t<<"\n";
                of<<pc->at(i).x<<" "<<pc->at(i).y<<"\n";
            }
            Point curr; 
            Point successor;
            of <<"e\n";
            for(size_t i = 0; i<mst.size(); i++) { 
                curr = pc->at(i);
                if(mst.at(i).successor) {
                    successor = pc->at(mst.at(i).successor->index);
                    //of<<curr.x<<" "<<curr.y<<" "<<curr.t<<"\n";
                    //of<<successor.x<<" "<<successor.y<<" "<<successor.t<<"\n\n\n";
                    of<<curr.x<<" "<<curr.y<<"\n";
                    of<<successor.x<<" "<<successor.y<<"\n\n\n";
                }
            }

            of <<"e\n";

            for(size_t i = 0; i<removed_edges.size(); i+=2) { 
                //of<<removed_edges.at(i).x<<" "<<removed_edges.at(i).y<<" "<<removed_edges.at(i).t<<"\n";
                //of<<removed_edges.at(i+1).x<<" "<<removed_edges.at(i+1).y<<" "<<removed_edges.at(i+1).t<<"\n\n\n";
                of<<removed_edges.at(i).x<<" "<<removed_edges.at(i).y<<"\n";
                of<<removed_edges.at(i+1).x<<" "<<removed_edges.at(i+1).y<<"\n\n\n";
            }

            of.close();
            return true;
        } else {
            of << "splot '-' using 1:2:3 with points palette pointtype 7 ps 1, '-' with lines lc 'black' title 'edges', '-' with lines lc 'green' title 'removed_edges'; pause -1;\n";
            for(size_t i = 0; i<pc->size(); i++) {
                of<<pc->at(i).x<<" "<<pc->at(i).y<<" "<<pc->at(i).t<<"\n";
            }
            Point curr;
            Point successor;
            of <<"e\n";
            for(size_t i = 0; i<mst.size(); i++) { 
                curr = pc->at(i);
                if(mst.at(i).successor) {
                    successor = pc->at(mst.at(i).successor->index);
                    of<<curr.x<<" "<<curr.y<<" "<<curr.t<<"\n";
                    of<<successor.x<<" "<<successor.y<<" "<<successor.t<<"\n\n\n";
                }
            }

            of <<"e\n";

            for(size_t i = 0; i<removed_edges.size(); i+=2) { 
                of<<removed_edges.at(i).x<<" "<<removed_edges.at(i).y<<" "<<removed_edges.at(i).t<<"\n";
                of<<removed_edges.at(i+1).x<<" "<<removed_edges.at(i+1).y<<" "<<removed_edges.at(i+1).t<<"\n\n\n";
            }

            of.close();
            return true;
        }   
    }






    UndirectedMST::UndirectedMST(PointCloud * pc, std::vector<int> * labels, double max_edge_dist, size_t max_branch_depth, bool edge_length_criterion, bool pca_angle_criterion, double pca_angle_value, int octree_target_size) {
        labels->clear();
        for(size_t i = 0; i<pc->size(); i++) {
            labels->push_back(-1);
            mst.push_back(UndirectedNode(i));
        }
        this->edge_length_criterion=edge_length_criterion;
        this->pca_angle_criterion=pca_angle_criterion;
        this->max_edge_dist=max_edge_dist;
        this->max_branch_depth=max_branch_depth;
        this->pca_angle_value=pca_angle_value;
        this->labels = labels;
        this->pc = pc;
        this->octree_target_size=octree_target_size;
        build_MST();
        //remove_time_gaps();
        //remove_long_edges();
        //remove_knees();

        //remove_branches();

        //color_nodes();
        //const char * fname = "pca_plot.gnuplot";
        //draw_pca(fname);
        //const char * fname2 = "plot_edges.gnuplot";
        //draw_pca_debug(fname2, false);
        //const char * fname3 = "pca_plot_edges2D.gnuplot";
        //draw_pca_debug(fname3, true);
    }

    // Custom comparison function 
    bool compareByDist(Edge a, Edge b) { 
        // Custom comparison logic 
        return a.dist < b.dist; // it sorts in ascending order 
    } 

    void UndirectedMST::build_MST() {
        size_t min_idx_i=std::numeric_limits<size_t>::max();
        size_t min_idx_j=std::numeric_limits<size_t>::max();
        float min_dist=std::numeric_limits<float>::max();
        float min_dist_small=std::numeric_limits<float>::max();
        size_t min_idx_small=std::numeric_limits<size_t>::max();
        size_t curr_idx=0;
        float curr_dist=0;
        size_t start_idx=0;
        std::vector<size_t> undiscovered;
        std::vector<ssize_t> nearest_neighbour_discovered;
        std::vector<size_t> discovered;

        const size_t n = mst.size();
        //float * dist_array = new float[mst.size()*mst.size()];
        std::vector<Edge> edges;
        edges.reserve((std::ceil((n*n)/2))-n);

        Edge temp;
        for(size_t i = 1; i<(n-1); i++) {
            for(size_t j = i+1; j<n; j++) {
                temp={i,j,pc->at(i).no_root_euclidian_distance(pc->at(j))};
                edges.push_back(temp);
            }
        }

        std::sort(edges.begin(), edges.end(), compareByDist);

        //std::unordered_map<size_t, 


        /* for(size_t i = 0; i<mst.size(); i++) {
            if(i!=start_idx) undiscovered.push_back(i);
        }
        discovered.push_back(start_idx);
        nearest_neighbour_discovered.push_back(-1);

        while(!undiscovered.empty()) {
            min_idx_i=std::numeric_limits<size_t>::max();
            min_idx_j=std::numeric_limits<size_t>::max();
            min_dist=std::numeric_limits<float>::max();
            for(size_t i=0; i<discovered.size(); i++) {
                if(nearest_neighbour_discovered.at(i)!=-1) {
                    curr_dist=dist_array[discovered.at(i)*n+nearest_neighbour_discovered.at(i)];
                    if(curr_dist<min_dist) {
                        min_dist=curr_dist;
                        min_idx_j=nearest_neighbour_discovered.at(i);
                        min_idx_i=i;
                    }
                } else {
                    //min_dist_small=std::numeric_limits<float>::max();
                    min_idx_small=std::numeric_limits<size_t>::max();
                    for(size_t j=0; j<undiscovered.size(); j++) {
                        curr_dist=dist_array[discovered.at(i)*n+undiscovered.at(j)];
                        if(curr_dist<min_dist) {
                            min_dist=curr_dist;
                            min_idx_j=j;
                            min_idx_i=i;
                        }
                        if(curr_dist<min_dist_small) {
                            //min_dist_small=curr_dist;
                            min_idx_small=j;
                        }
                    }
                    nearest_neighbour_discovered.at(i)=undiscovered.at(min_idx_small);
                }
            }
            mst.at(discovered.at(min_idx_i)).neighbours.push_back(&mst.at(undiscovered.at(min_idx_j)));
            mst.at(undiscovered.at(min_idx_j)).neighbours.push_back(&mst.at(discovered.at(min_idx_i)));
            discovered.push_back(undiscovered.at(min_idx_j));
            nearest_neighbour_discovered.push_back(-1);
            for(size_t i = 0; i<nearest_neighbour_discovered.size(); i++) {
                if(nearest_neighbour_discovered.at(i)==undiscovered.at(min_idx_j)) nearest_neighbour_discovered.at(i)=-1;
            }
            undiscovered.erase(undiscovered.begin() + min_idx_j);
        } */



        //delete [] dist_array;

    }



    bool UndirectedMST::draw_mst(const char * fname) {
        // saves cloud in gnuplot file
        std::ofstream of(fname);
        of << std::fixed;  // set float style
        if (!of.is_open()) {
        std::cerr << "[Error] could not save under '" << fname << std::endl;
        return false;
        }

        of << "splot '-' using 1:2:3 with points lc 'black' pointtype 7 ps 1, '-' with lines palette title 'edges'; pause -1;\n";
        for(size_t i = 0; i<pc->size(); i++) {
            of<<pc->at(i).x<<" "<<pc->at(i).y<<" "<<pc->at(i).t<<"\n";
        }
        Point curr;
        Point successor;
        of <<"e\n";
        for(size_t i = 0; i<mst.size(); i++) { 
            curr = pc->at(i);
            for(size_t j = 0; i<mst.at(i).neighbours.size(); i++) {
                successor = pc->at(mst.at(i).neighbours.at(j)->index);
                of<<curr.x<<" "<<curr.y<<" "<<curr.t<<"\n";
                of<<successor.x<<" "<<successor.y<<" "<<successor.t<<"\n\n\n";
            }
        }
        of <<"e\n";
        of.close();
        return true;
    }


    