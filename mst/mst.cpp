    #include "mst.hpp"
    #include <sstream>
    #include <fstream>
    #include <limits>
    #include <deque>
    #include <unordered_map>
    #include <utility>
    #include "../utils/mathstuff.hpp"
    #include "emst.hpp"
    #include <iostream>


    UndirectedMST::UndirectedMST(PointCloud * pc, std::vector<int> * labels, float max_edge_dist, size_t min_cluster_size, const size_t ndims, double max_mean_edge_len) {
        this->max_mean_edge_len = max_mean_edge_len;
        this->min_cluster_size = min_cluster_size;
        this->max_edge_dist = max_edge_dist;
        this->pc = pc;
        this->labels = labels;
        this->n = pc->size();
        labels->clear();
        for(size_t i = 0; i<n; i++) {
            labels->push_back(-1);
            mst.push_back(UndirectedNode(i));
        }

        // calc MST with the library Euclidian MST by AndrewB330
        Point0 curr;
        if(ndims==3) {
            std::vector<Point<3>> points(n);
            for(size_t i = 0; i<n; i++) {
                curr = pc->at(i);
                for(size_t k = 0; k<(ndims-1); k++) {
                    points[i][k] = curr.data.at(k);
                }
                points[i][ndims-1] = curr.t;
            }
            KdTreeSolver<3> solver(points);
            double total_length = solver.get_total_length();
            std::vector<Edge> edges = solver.get_solution();

            // turn Edges into my own MST nodes because further processing easier that way  
            Point0 a=Point0(ndims);
            Point0 b=Point0(ndims);
            float dist;
            //double mean=0.0;
            for(size_t i = 0; i<edges.size(); i++) {
                a = pc->at(std::get<0>(edges.at(i)));
                b = pc->at(std::get<1>(edges.at(i)));
                dist = a.no_root_euclidian_distance(b);
                //mean += dist;
                mst.at(std::get<0>(edges.at(i))).neighbours.push_back(&mst.at(std::get<1>(edges.at(i))));
                mst.at(std::get<0>(edges.at(i))).distances.push_back(dist);
                mst.at(std::get<1>(edges.at(i))).neighbours.push_back(&mst.at(std::get<0>(edges.at(i))));
                mst.at(std::get<1>(edges.at(i))).distances.push_back(dist);
            }
            //mean/=edges.size();
            //double sd=0.0;
            //for(size_t i = 0; i<edges.size(); i++) {
            //    a = pc->at(std::get<0>(edges.at(i)));
            //    b = pc->at(std::get<1>(edges.at(i)));
            //    dist = a.no_root_euclidian_distance(b);
            //    sd += (dist-mean)*(dist-mean);
            //}
            //sd = std::sqrt(sd/(edges.size()-1));
            //std::cout<<mean<<","<<sd<<","<<pc->size()<<"\n";
            remove_long_edges(max_edge_dist);
            //remove_local_long_edges(20, 1.0);
        } else if(ndims==4) {
            std::vector<Point<4>> points(n);
            for(size_t i = 0; i<n; i++) {
                curr = pc->at(i);
                for(size_t k = 0; k<(ndims-1); k++) {
                    points[i][k] = curr.data.at(k);
                }
                points[i][ndims-1] = curr.t;
            }
            KdTreeSolver<4> solver(points);
            double total_length = solver.get_total_length();
            std::vector<Edge> edges = solver.get_solution();

            // turn Edges into my own MST nodes because further processing easier that way  
            Point0 a=Point0(ndims);
            Point0 b=Point0(ndims);
            float dist;
            for(size_t i = 0; i<edges.size(); i++) {
                a = pc->at(std::get<0>(edges.at(i)));
                b = pc->at(std::get<1>(edges.at(i)));
                dist = a.no_root_euclidian_distance(b);
                mst.at(std::get<0>(edges.at(i))).neighbours.push_back(&mst.at(std::get<1>(edges.at(i))));
                mst.at(std::get<0>(edges.at(i))).distances.push_back(dist);
                mst.at(std::get<1>(edges.at(i))).neighbours.push_back(&mst.at(std::get<0>(edges.at(i))));
                mst.at(std::get<1>(edges.at(i))).distances.push_back(dist);
            }
            remove_long_edges(max_edge_dist);
        } else {
            std::cout<<"[Error in UndirectedMST::UndirectedMST()]: wrong number of dimensions, only 3 or 4 allowed\n";
        }
    }


    UndirectedMST::UndirectedMST(PointCloud* pc, std::vector<int> * labels, size_t min_depth, size_t min_cluster_size, const size_t ndims) {
        this->min_cluster_size = min_cluster_size;
        this->num_clusters = -1;
        this->min_depth = min_depth;
        this->pc = pc;
        this->labels = labels;
        this->n = pc->size();
        labels->clear();
        for(size_t i = 0; i<n; i++) {
            labels->push_back(-1);
            mst.push_back(UndirectedNode(i));
        }

        // calc MST with the library Euclidian MST by AndrewB330
        if(ndims==3) {
            std::vector<Point<3>> points(n);
            Point0 curr;
            for(size_t i = 0; i<n; i++) {
                curr = pc->at(i);
                //points[i][0] = curr.x;
                //points[i][1] = curr.y;
                //points[i][2] = curr.t;
                for(size_t k = 0; k<(ndims-1); k++) {
                    points[i][k] = curr.data.at(k);
                }
                points[i][ndims-1] = curr.t;
            }
            KdTreeSolver<3> solver(points);
            double total_length = solver.get_total_length();
            std::vector<Edge> edges = solver.get_solution();
           
            // turn Edges into my own MST nodes because further processing easier that way  
            Point0 a;
            Point0 b;
            float dist;
            for(size_t i = 0; i<edges.size(); i++) {
                //std::cout<<std::get<0>(edges.at(i))<<" "<<std::get<1>(edges.at(i))<<"\n";
                a = pc->at(std::get<0>(edges.at(i)));
                b = pc->at(std::get<1>(edges.at(i)));
                dist = a.no_root_euclidian_distance(b);
                mst.at(std::get<0>(edges.at(i))).neighbours.push_back(&mst.at(std::get<1>(edges.at(i))));
                mst.at(std::get<0>(edges.at(i))).distances.push_back(dist);
                mst.at(std::get<1>(edges.at(i))).neighbours.push_back(&mst.at(std::get<0>(edges.at(i))));
                mst.at(std::get<1>(edges.at(i))).distances.push_back(dist);
            }
            remove_branches();
        } else if(ndims==4) {
            std::vector<Point<4>> points(n);
            Point0 curr;
            for(size_t i = 0; i<n; i++) {
                curr = pc->at(i);
                //points[i][0] = curr.x;
                //points[i][1] = curr.y;
                //points[i][2] = curr.t;
                for(size_t k = 0; k<(ndims-1); k++) {
                    points[i][k] = curr.data.at(k);
                }
                points[i][ndims-1] = curr.t;
            }
            KdTreeSolver<4> solver(points);
            double total_length = solver.get_total_length();
            std::vector<Edge> edges = solver.get_solution();
           
            // turn Edges into my own MST nodes because further processing easier that way  
            Point0 a;
            Point0 b;
            float dist;
            for(size_t i = 0; i<edges.size(); i++) {
                //std::cout<<std::get<0>(edges.at(i))<<" "<<std::get<1>(edges.at(i))<<"\n";
                a = pc->at(std::get<0>(edges.at(i)));
                b = pc->at(std::get<1>(edges.at(i)));
                dist = a.no_root_euclidian_distance(b);
                mst.at(std::get<0>(edges.at(i))).neighbours.push_back(&mst.at(std::get<1>(edges.at(i))));
                mst.at(std::get<0>(edges.at(i))).distances.push_back(dist);
                mst.at(std::get<1>(edges.at(i))).neighbours.push_back(&mst.at(std::get<0>(edges.at(i))));
                mst.at(std::get<1>(edges.at(i))).distances.push_back(dist);
            }
            remove_branches();
        } else {
            std::cout<<"[Error in UndirectedMST::UndirectedMST()]: wrong number of dimensions, only 3 or 4 allowed\n";
        }
    }


    size_t UndirectedMST::colorMST() {
        return color_nodes();
    }


    bool UndirectedMST::draw_mst(const char * fname, bool is4d) {
        std::ofstream of(fname);
        of << std::fixed;  // set float style
        if (!of.is_open()) {
            std::cerr << "[Error] could not save under '" << fname << std::endl;
            return false;
        }

        //of << "unset colorbox\nunset border\nunset xtics\nunset ytics\nunset ztics\nset style arrow 1 head filled size screen 0.02,15,30 lw 2 lc 'black'\n";

        /* double min_dist=std::numeric_limits<double>::max();
        size_t min_idx=0;
        for(size_t i = 0; i<n; i++) {
            of << "set arrow from ";
            for(size_t m = 0; m<pc->at(i).data.size(); m++) {
                of<<pc->at(i).data.at(m)<<",";
            }
            of << pc->at(i).t << " to ";
            min_dist=std::numeric_limits<double>::max();
            size_t j=i-1;
            while((j<std::numeric_limits<size_t>::max()) && (pc->at(j).t >= pc->at(i).t)) {
                if(pc->at(i).no_root_euclidian_distance(pc->at(j)) < min_dist) {
                    min_dist = pc->at(i).no_root_euclidian_distance(pc->at(j));
                    min_idx = j;
                }
                j--;
            }
            for(j=i+1; j<n; j++) {
                if(pc->at(i).no_root_euclidian_distance(pc->at(j)) < min_dist) {
                    min_dist = pc->at(i).no_root_euclidian_distance(pc->at(j));
                    min_idx = j;
                }
            }
            for(size_t m = 0; m<pc->at(min_idx).data.size(); m++) {
                of<<pc->at(min_idx).data.at(m)<<",";
            }
            of << pc->at(min_idx).t << " as 1\n";
        } */

        //of << "splot '-' using 1:2:3 with points lc 'black' pointtype 7 ps 1, '-' with lines palette title 'edges'; pause -1;\n";
        of << "unset colorbox\nunset border\nunset xtics\nunset ytics\nunset ztics\nsplot '-' using 1:2:3 with points palette pointtype 7 ps 2, '-' with lines lc 'black' title 'edges'; pause -1;\n";
        //of << "splot '-' using 1:2:3 with points palette pointtype 7 ps 1; pause -1;\n";
        if(is4d) {
            for(size_t i = 0; i<pc->size(); i++) {
                for(size_t k = 0; k<pc->at(i).data.size(); k++) {
                    of<<pc->at(i).data.at(k)<<" ";
                }
                of << "\n";
            }
        } else {
            for(size_t i = 0; i<pc->size(); i++) {
                for(size_t k = 0; k<pc->at(i).data.size(); k++) {
                    of<<pc->at(i).data.at(k)<<" ";
                }
                of << pc->at(i).t<<"\n";
            }
        }
        size_t k;
        size_t l;
        of <<"e\n";

        for(size_t i = 0; i<n; i++) { 
            for(size_t j = 0; j<mst.at(i).neighbours.size(); j++) {
                k=mst.at(i).index;
                l=mst.at(i).neighbours.at(j)->index;
                if(is4d) {
                    for(size_t m = 0; m<pc->at(i).data.size(); m++) {
                        of<<pc->at(k).data.at(m)<<" ";
                    }
                    of << "\n";
                    for(size_t m = 0; m<pc->at(i).data.size(); m++) {
                        of<<pc->at(l).data.at(m)<<" ";
                    }
                    of << "\n\n\n";
                } else {
                    for(size_t m = 0; m<pc->at(i).data.size(); m++) {
                        of<<pc->at(k).data.at(m)<<" ";
                    }
                    of << pc->at(k).t << "\n";
                    for(size_t m = 0; m<pc->at(i).data.size(); m++) {
                        of<<pc->at(l).data.at(m)<<" ";
                    }
                    of << pc->at(l).t << "\n\n\n";
                }
            }
        }

        of.close();
        return true;
    }

    float UndirectedMST::get_edge_len(std::pair<size_t, size_t> edge) {
        Point0 a = pc->at(std::get<0>(edge));
        Point0 b = pc->at(std::get<1>(edge));
        return a.no_root_euclidian_distance(b);
    }


    void UndirectedMST::remove_long_edges(double max_edge_len) {
        std::vector<UndirectedNode*> out_neighbours;
        std::vector<float> out_dist;
        for(size_t i = 0; i<n; i++) {
            out_neighbours.clear();
            out_dist.clear();
            for(size_t j = 0; j<mst.at(i).neighbours.size(); j++) {
                if(mst.at(i).distances.at(j) <= max_edge_len) {
                    out_neighbours.push_back(mst.at(i).neighbours.at(j));
                    out_dist.push_back(mst.at(i).distances.at(j));
                }
            }
            mst.at(i).distances=out_dist;
            mst.at(i).neighbours=out_neighbours;
        }
    }

    size_t UndirectedMST::color_nodes() {
        size_t n = pc->size();
        int cluster_id = 0;
        double mean_edge_dist=0;
        std::vector<size_t> count_label_members;
        std::vector<double> mean_edge_dists;
        //double max_mean_edge_dist=100.0;
        // initial labelling
        for(size_t i = 0; i<n; i++) {
            if(labels->at(i)==-1) {
                mean_edge_dist = expand(&mst.at(i), cluster_id);
                mean_edge_dists.push_back(mean_edge_dist);
                cluster_id++;
                //if(cluster_id==collision) cluster_id++;
                //else merge_labels(collision, cluster_id);
            }
        }

        double count[cluster_id]={0}; //assign all clusters with small amount of data the label noise (-1)
        for(size_t i = 0; i<n; i++) {
            count[labels->at(i)]++;
        }
        if(max_mean_edge_len>0) {
            for(size_t i = 0; i<cluster_id; i++) {
                if((count[i]<min_cluster_size) || (mean_edge_dists[i] > max_mean_edge_len)) {
                    for(size_t j = 0; j<n; j++) {
                        if(labels->at(j)==i) labels->at(j) = -1;
                    }
                    count[i]=0;
                }
            }
        } else {
            for(size_t i = 0; i<cluster_id; i++) {
                if((count[i]<min_cluster_size)) {
                    for(size_t j = 0; j<n; j++) {
                        if(labels->at(j)==i) labels->at(j) = -1;
                    }
                    count[i]=0;
                }
            }
        }
        // some labels are not consecutive anymore
        size_t number_removed_clusters=0;
        for(size_t i = 0; i<cluster_id; i++){
            if(count[i]==0) number_removed_clusters++;
            else {
                for(size_t j = 0; j<n; j++) {
                    if(labels->at(j)==i) labels->at(j) -= number_removed_clusters;
                }
            }
        }
        num_clusters = cluster_id;

        return cluster_id-number_removed_clusters; 
    }


    /* void UndirectedMST::merge_labels(size_t label_a, size_t label_b) {
        size_t label_to_overwrite = (label_a > label_b) ? label_a : label_b;
        size_t label_left_standing = (label_a > label_b) ? label_b : label_a;
        for(size_t i = 0; i<n; i++) {
            if(labels->at(i) == label_to_overwrite) labels->at(i) = label_left_standing;
        }
    } */


    double UndirectedMST::expand(UndirectedNode * node, int cluster_id) {
        std::deque<UndirectedNode*> open;
        std::unordered_map<size_t, bool> closed;
        //std::vector<double> dists;
        open.push_front(node);
        UndirectedNode* curr_node;
        double mean_edge_dist=0.0;
        size_t counter=0;
        while(!open.empty()) {
            curr_node = open.front();
            open.pop_front();
            closed[curr_node->index] = true;
            labels->at(curr_node->index) = cluster_id;
            for(size_t i = 0; i<curr_node->neighbours.size(); i++) {
                if(!closed[curr_node->neighbours.at(i)->index]) {
                    //if(labels->at(curr_node->neighbours.at(i)->index) != -1) return labels->at(curr_node->neighbours.at(i)->index);
                    open.push_front(curr_node->neighbours.at(i));
                    mean_edge_dist+=curr_node->distances.at(i);
                    //dists.push_back(curr_node->distances.at(i));
                    counter++;
                }
            }
        }
        mean_edge_dist /= counter;

        return mean_edge_dist;
    }


    void UndirectedMST::edge_dist_at_depth(UndirectedNode * node, UndirectedNode * taboo_node, size_t depth, double* mean, double* sd) {
        std::deque<UndirectedNode*> open;
        std::deque<size_t> depth_mem;
        std::unordered_map<size_t, bool> closed;
        std::vector<double> dists;
        open.push_front(node);
        depth_mem.push_front(0);
        UndirectedNode* curr_node;
        closed[taboo_node->index] = true;  // don't allow taboo node
        double mean_edge_dist=0.0;
        size_t counter=0;
        size_t curr_depth=0;
        while(!open.empty()) {
            curr_node = open.front();
            open.pop_front();
            curr_depth = depth_mem.front();
            depth_mem.pop_front();
            closed[curr_node->index] = true;
            if(curr_depth<depth) {
                for(size_t i = 0; i<curr_node->neighbours.size(); i++) {
                    if(!closed[curr_node->neighbours.at(i)->index]) {
                        //if(labels->at(curr_node->neighbours.at(i)->index) != -1) return labels->at(curr_node->neighbours.at(i)->index);
                        open.push_front(curr_node->neighbours.at(i));
                        depth_mem.push_front(curr_depth+1);
                        mean_edge_dist+=curr_node->distances.at(i);
                        counter++;
                        dists.push_back(curr_node->distances.at(i));
                    }
                }
            }
        }
        mean_edge_dist /= counter;
        *mean = mean_edge_dist;
        double std=0.0;
        double curr=0.0;
        for(size_t i = 0; i<dists.size(); i++) {
            curr = (dists.at(i)-mean_edge_dist);
            std += curr*curr;
        }
        *sd = std::sqrt(std/(dists.size()-1));
    }

    void UndirectedMST::remove_local_long_edges(size_t depth, double c) {
        double mean1=0.0;
        double sd1=0.0;
        double mean2=0.0;
        double sd2=0.0;

        double threshold=0.0;

        std::vector<UndirectedNode*> out_neighbours;
        std::vector<float> out_dist;
        for(size_t i = 0; i<n; i++) {
            out_neighbours.clear();
            out_dist.clear();
            for(size_t j = 0; j<mst.at(i).neighbours.size(); j++) {
                edge_dist_at_depth(&mst.at(i), mst.at(i).neighbours.at(j), depth, &mean1, &sd1);
                edge_dist_at_depth(mst.at(i).neighbours.at(j), &mst.at(i), depth, &mean2, &sd2);
                threshold = std::max((mean1+c*sd1), (mean2+c*sd2));
                if(!std::isfinite(mean1+c*sd1)) threshold = (mean2+c*sd2);
                if(!std::isfinite(mean2+c*sd2)) threshold = (mean1+c*sd1);
                if(!std::isfinite(threshold)) threshold = max_edge_dist;

                //if((mean1+c*sd1)>(mean2+c*sd2)) std::cout<<(mst.at(i).distances.at(j)-mean1)/sd1<<"\n";
                //else std::cout<<(mst.at(i).distances.at(j)-mean2)/sd2<<"\n";

                if(mst.at(i).distances.at(j) <= threshold) {
                    out_neighbours.push_back(mst.at(i).neighbours.at(j));
                    out_dist.push_back(mst.at(i).distances.at(j));
                }
            }
            mst.at(i).distances=out_dist;
            mst.at(i).neighbours=out_neighbours;
        }

    }

    bool UndirectedMST::probe_depth(size_t node_index, size_t node_neighbour, size_t max_depth) {
        UndirectedNode node = mst.at(node_neighbour);
        std::deque<UndirectedNode*> open;
        std::deque<size_t> depth;
        std::unordered_map<size_t, bool> closed;
        closed[node_index] = true;
        open.push_front(&node);
        depth.push_back(1);
        UndirectedNode* curr_node;
        size_t curr_depth;
        float min_dist=std::numeric_limits<float>::max();
        ssize_t min_idx=-1;
        while(!open.empty()) {
            curr_node = open.front();
            curr_depth = depth.front();
            //std::cout<<curr_node->index<<" "<<curr_depth<<"\n";
            open.pop_front();
            depth.pop_front();
            if(curr_depth>max_depth) return true;
            closed[curr_node->index] = true;
            //min_dist=std::numeric_limits<float>::max();
            //min_idx=-1;
            for(size_t i = 0; i<curr_node->neighbours.size(); i++) {
                if(!closed[curr_node->neighbours.at(i)->index]) {
                    open.push_front(curr_node->neighbours.at(i));
                    depth.push_back(curr_depth+1);
                }
                /* if(curr_node->distances.at(i) < min_dist) {
                    min_dist = curr_node->distances.at(i);
                    min_idx = i;
                } */
            }
            /* if(min_idx!=-1) {
                open.push_front(curr_node->neighbours.at(min_idx));
                depth.push_back(curr_depth+1);
            } */
        }
        return false;
    }


    void UndirectedMST::remove_branches() {
        size_t num_branches=0;
        std::vector<UndirectedNode*> new_neighbour_list;
        std::vector<float> new_neighbour_dists;
        float shortest_branch=std::numeric_limits<float>::max();
        size_t shortest_branch_idx=0;
        std::vector<std::pair<size_t, size_t>> to_remove; // first element is index of node where the edges are cut, second one the index of the node that the shortest edge goes to 
        std::vector<float> dist_to_stay;  // distance for every pair in to_remove... 

        size_t i;

        for(i = 0; i<mst.size(); i++) { 
            num_branches=0;
            if(mst.at(i).neighbours.size()>2) {
                shortest_branch=std::numeric_limits<float>::max();
                for(size_t j = 0; j<mst.at(i).neighbours.size(); j++) {
                    if(probe_depth(i,mst.at(i).neighbours.at(j)->index, min_depth)) {
                       num_branches++;
                    } 
                    if(shortest_branch>mst.at(i).distances.at(j)) {
                        shortest_branch=mst.at(i).distances.at(j);
                        shortest_branch_idx=mst.at(i).neighbours.at(j)->index;
                    }
                }
                if(num_branches>2) {

                    to_remove.push_back(std::make_pair(i, shortest_branch_idx));
                    dist_to_stay.push_back(shortest_branch);

                    //for(size_t j = 0; j<mst.at(i).neighbours.size(); j++) {
                    //    new_neighbour_list.clear();
                    //    new_neighbour_dists.clear();
                    //    for(size_t k = 0; k<mst.at(i).neighbours.at(j)->neighbours.size(); k++) {
                    //        if(mst.at(i).neighbours.at(j)->neighbours.at(k)->index!=i) {
                    //            new_neighbour_list.push_back(mst.at(i).neighbours.at(j)->neighbours.at(k));
                    //            new_neighbour_dists.push_back(mst.at(i).neighbours.at(j)->distances.at(k));
                    //        }
                    //    }
                    //    mst.at(i).neighbours.at(j)->neighbours = new_neighbour_list;
                    //    mst.at(i).neighbours.at(j)->distances = new_neighbour_dists;
                    //}
                    //mst.at(i).neighbours.clear();
                    //mst.at(i).distances.clear();
                    //mst.at(i).neighbours.push_back(&mst.at(shortest_branch_idx));
                    //mst.at(i).distances.push_back(shortest_branch);
                    //mst.at(shortest_branch_idx).neighbours.push_back(&mst.at(i));
                    //mst.at(shortest_branch_idx).distances.push_back(shortest_branch);
                }   
            }
        }

        for(size_t l = 0; l<to_remove.size(); l++) {
            i = to_remove.at(l).first;
            shortest_branch_idx = to_remove.at(l).second;
            shortest_branch = dist_to_stay.at(l);
            for(size_t j = 0; j<mst.at(i).neighbours.size(); j++) {
                    new_neighbour_list.clear();
                    new_neighbour_dists.clear();
                    for(size_t k = 0; k<mst.at(i).neighbours.at(j)->neighbours.size(); k++) {
                        if(mst.at(i).neighbours.at(j)->neighbours.at(k)->index!=i) {
                            new_neighbour_list.push_back(mst.at(i).neighbours.at(j)->neighbours.at(k));
                            new_neighbour_dists.push_back(mst.at(i).neighbours.at(j)->distances.at(k));
                        }
                    }
                    mst.at(i).neighbours.at(j)->neighbours = new_neighbour_list;
                    mst.at(i).neighbours.at(j)->distances = new_neighbour_dists;
            }
            mst.at(i).neighbours.clear();
            mst.at(i).distances.clear();
            mst.at(i).neighbours.push_back(&mst.at(shortest_branch_idx));
            mst.at(i).distances.push_back(shortest_branch);
            mst.at(shortest_branch_idx).neighbours.push_back(&mst.at(i));
            mst.at(shortest_branch_idx).distances.push_back(shortest_branch);
        }


    }
    

    /* PointCloud n1;
    PointCloud n2;
    // für test_kreuzung2 mit octree target_size 10
    //int idx123=10;     // an dem Knick
    //int idx123=85;     // an einer anderen inkonsistenten Kante
    //int idx123=76;
    //int idx123=30;       // innerhalb einer langen Strecke 
    //int idx123=98;      // bei big boy small mit octree target size 10 eine inkonsistente Kante und ein Beispiel für unrobustheit von PCA 
    //int idx123=393; //390   // bei input_test_small ein Beispiel für n1_problem
    //int idx123=10;
    //int idx123=144;        // innerhalb einer gebogenen Strecke
    //int idx123=101;
    //int idx123=145;
    //int idx123=55;
    //int idx123=3802;
    int idx123=48;
    PointCloud pca_idx123_n1;
    PointCloud pca_idx123_n2; 


    void RANSAC(std::vector<Point0> * neighbourhood, float t, size_t * i, size_t * j) {
        MathStuff ms;
        size_t max_i=0;
        size_t max_j=0;
        size_t max_num_supporters=0;
        size_t curr_supporters=0;
        Point0 a;
        Point0 b;
        Point0 c;
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
                if(curr_supporters>=max_num_supporters) {
                    max_i=i;
                    max_j=j;
                    max_num_supporters=curr_supporters;
                }
            }
        }
        *i=max_i;
        *j=max_j;
    }*/


    DirectedMST::DirectedMST(PointCloud * pc, std::vector<int> * labels, double max_edge_dist, int octree_target_size, size_t min_cluster_size) {

        labels->clear();
        for(size_t i = 0; i<pc->size(); i++) {
            labels->push_back(-1);
            mst.push_back(DirectedNode(i));
        }
        this->max_edge_dist = max_edge_dist;
        this->labels = labels;
        this->pc = pc;
        this->octree_target_size = octree_target_size;
        this->min_cluster_size = min_cluster_size;
        build_MST();
        remove_long_edges();

        //color_nodes();
    }

    void DirectedMST::build_MST() { // NO FLANN
        size_t n = pc->size();
        Point0 curr;
        const size_t big_size_t=std::numeric_limits<size_t>::max();
        const double big_double=std::numeric_limits<double>::max();
        double min_distance=big_double;
        size_t min_index=big_size_t;
        double curr_distance=0;
        double delta_t_squared=0;
        size_t j;

        for(size_t i = 0; i<n; i++) {
            curr = pc->at(i);
            min_distance=big_double;
            min_index=big_size_t;


            /* if(((i+1) < pc->size()) && !(pc->at(i+1).t != pc->at(i).t)) {
                j=i-1;
                while((j<std::numeric_limits<size_t>::max()) && (pc->at(j).t >= pc->at(i).t)) {
                    if(pc->at(i).no_root_euclidian_distance(pc->at(j)) < min_distance) {
                        min_distance = pc->at(i).no_root_euclidian_distance(pc->at(j));
                        min_index = j;
                    }
                    j--;
                }
            } */

            for(j = i+1; j<n; j++) {
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




    void DirectedMST::remove_long_edges() {
        //std::vector<size_t> snipsnap;  
        size_t a, b;
        //double global_max_dist_threshold = choose_global_max_dist_threshold();
        double global_max_dist_threshold = max_edge_dist;
        double local_max_dist_threshold=global_max_dist_threshold;
        //std::cout<<"Chose global max_dist_threshold="<<std::sqrt(global_max_dist_threshold)<<"\n";
        for(size_t i = 0; i<mst.size(); i++) {
            //if(this->max_edge_dist==-2) local_max_dist_threshold=choose_local_max_dist_threshold(i);
            if(mst.at(i).successor) {
                //std::cout<<mst.at(i).dist_successor<<std::endl;
                //std::cout<<max_edge_dist<<"="<<global_max_dist_threshold<<"="<<local_max_dist_threshold<<std::endl;
                if((mst.at(i).dist_successor > local_max_dist_threshold)) {
                    delete_edge(i);
                }
            }
        }
    }

    /*
    void DirectedMST::remove_knees() {
        const float ransac_t=200.0;
        const size_t max_depth=5;
        //size_t curr_depth=0;
        PointCloud vectors_n1;
        PointCloud vectors_n2;
        std::vector<Point0> ps;
        Point0 a_n1;
        Point0 b_n1;
        Point0 a_n2;
        Point0 b_n2;
        MathStuff mathstuff;
        double cos_similarity_n1=0.0;
        double cos_similarity_n2=0.0;
        for(size_t index = 0; index<mst.size(); index++) {
            if(mst.at(index).successor && !((pc->at(index)-pc->at(mst.at(index).successor->index)) == Point0(0,0,0))) { //&& isBranch(index)) {

                //std::cout<<index<<std::endl;

                BFS(index, &vectors_n1, &vectors_n2, max_depth);

                if((vectors_n1.size()>=(size_t)(max_depth/2)) && (vectors_n2.size()>=(size_t)(max_depth/2))) {
                    //mathstuff.orthogonal_lsq(vectors_n1, &a_n1, &b_n1);
                    //mathstuff.orthogonal_lsq(vectors_n2, &a_n2, &b_n2);
                    size_t idx_i;
                    size_t idx_j;
                    RANSAC(&vectors_n1, ransac_t, &idx_i, &idx_j);
                    a_n1 = pc->at(index);
                    b_n1 = vectors_n1.at(idx_i)-vectors_n1.at(idx_j);
                    RANSAC(&vectors_n2, ransac_t, &idx_i, &idx_j);
                    a_n2 = pc->at(mst.at(index).successor->index);
                    b_n2 = vectors_n2.at(idx_i)-vectors_n2.at(idx_j);
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
                    }
                } else if(vectors_n1.size()>=(size_t)(max_depth/2)) {
                    //mathstuff.orthogonal_lsq(vectors_n1, &a_n1, &b_n1);
                    size_t idx_i;
                    size_t idx_j;
                    RANSAC(&vectors_n1, ransac_t, &idx_i, &idx_j);
                    a_n1 = pc->at(index);
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
                    } 
                } else if(vectors_n2.size()>=(size_t)(max_depth/2)) {
                    //mathstuff.orthogonal_lsq(vectors_n2, &a_n2, &b_n2);
                    size_t idx_i;
                    size_t idx_j;
                    RANSAC(&vectors_n2, ransac_t, &idx_i, &idx_j);
                    a_n2 = pc->at(mst.at(index).successor->index);
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
                    }
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
    } */      


    ssize_t DirectedMST::BFS(size_t index, PointCloud * n1, PointCloud * n2, size_t max_depth) {
        n1->clear();
        n2->clear();
        size_t curr_depth=0;
        DirectedNode * curr_node;
        std::deque<DirectedNode*> open;
        std::deque<size_t> depth;
        std::unordered_map<size_t, bool> closed;

        DirectedNode* min_node;
        double min_distance=std::numeric_limits<double>::max();

        // !!!!!!!!!!!!!!!!!! n1 !!!!!!!!!!!!!!!!!!!
        for(size_t i = 0; i<mst.at(index).predecessor.size(); i++) {
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
                min_distance=std::numeric_limits<double>::max();
                for(size_t i = 0; i<curr_node->predecessor.size(); i++) {
                    if(curr_node->dist_predecessor.at(i) < min_distance) {
                        min_node = curr_node->predecessor.at(i);
                        min_distance=curr_node->dist_predecessor.at(i);
                    }
                    //open.push_back(curr_node->predecessor.at(i));
                    //depth.push_back(curr_depth+1);
                }
                open.push_back(min_node);
                depth.push_back(curr_depth+1);
            }
        }

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
            }
            closed[curr_node->index] = true;
        }

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
            double val=6*std::sqrt(2*(double)(octree_target_size*octree_target_size));
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

    /*
    bool DirectedMST::isBranch(size_t i) {
        bool remove_edge=false;
        size_t n = mst.size();
        if(mst.at(i).predecessor.size() >= 2) {
            for(size_t j = 0; j < mst.at(i).predecessor.size(); j++) {
                remove_edge=DFS(mst.at(i).predecessor.at(j),0);
                if(!remove_edge) {
                    return false;
                }
            }
        } else return false;   
        return true;
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
                if(to_be_removed.size()==mst.at(i).predecessor.size()) {
                for(size_t j : to_be_removed) {
                    mst.at(i).predecessor.at(j)->successor = NULL;
                }
                }
            }
        } 
    }*/

    /* bool DirectedMST::DFS(DirectedNode * n, size_t curr_depth) {
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
    } */


    int DirectedMST::expand(DirectedNode * node, int cluster_id) {
        std::deque<DirectedNode*> open;
        std::unordered_map<size_t, bool> closed;
        open.push_front(node);
        DirectedNode* curr_node;
        while(!open.empty()) {
            curr_node = open.front();
            open.pop_front();
            closed[curr_node->index] = true;
            labels->at(curr_node->index) = cluster_id;
            for(size_t i = 0; i<curr_node->predecessor.size(); i++) {
                if(!closed[curr_node->predecessor.at(i)->index]) {
                    //if(labels->at(curr_node->neighbours.at(i)->index) != -1) return labels->at(curr_node->neighbours.at(i)->index);
                    open.push_front(curr_node->predecessor.at(i));
                }
            }
            if((curr_node->successor) && (!closed[curr_node->successor->index])) open.push_front(curr_node->successor);
        }

        return cluster_id;
        /* if(labels->at(node->index) != -1) {
            return labels->at(node->index);
        }
        int collision=cluster_id;
        if(node->successor) 
            collision = expand(node->successor, cluster_id);  
        labels->at(node->index) = collision;
        return collision; */
    }


    size_t DirectedMST::colorMST() {
        size_t n = pc->size();
        int cluster_id = 0;
        int collision=0;
        for(size_t i = 0; i<n; i++) {
            if(labels->at(i)==-1) {
                collision = expand(&mst.at(i), cluster_id);
                cluster_id++; 
            }
        }

        ////// labels with small amounts of data are noise?

        double count[cluster_id]={0}; //assign all clusters with small amount of data the label noise (-1)
        for(size_t i = 0; i<n; i++) {
            count[labels->at(i)]++;
        }
        for(size_t i = 0; i<cluster_id; i++) {
            if((count[i]<min_cluster_size)) {
                for(size_t j = 0; j<n; j++) {
                    if(labels->at(j)==i) labels->at(j) = -1;
                }
                count[i]=0;
            }
        }
        // some labels are not consecutive anymore
        size_t number_removed_clusters=0;
        for(size_t i = 0; i<cluster_id; i++){
            if(count[i]==0) number_removed_clusters++;
            else {
                for(size_t j = 0; j<n; j++) {
                    if(labels->at(j)==i) labels->at(j) -= number_removed_clusters;
                }
            }
        }
        size_t num_clusters = cluster_id;

        return cluster_id-number_removed_clusters; 


    }

    
    bool DirectedMST::draw_mst(const char * fname, bool is4d) {
        size_t version=0;
        // saves cloud in gnuplot file
        std::ofstream of(fname);
        of << std::fixed;  // set float style
        if (!of.is_open()) {
            std::cerr << "[Error] could not save under '" << fname << std::endl;
            return false;
        }

        of << "unset colorbox\nunset border\nunset xtics\nunset ytics\nunset ztics\nset style arrow 1 head filled size screen 0.02,15,30 lw 2 lc 'black'\n";


        Point0 curr;
        Point0 successor;
        for(size_t i = 0; i<mst.size(); i++) { 
            curr = pc->at(i);
            if(mst.at(i).successor) {
                of << "set arrow from ";
                successor = pc->at(mst.at(i).successor->index);
                for(size_t m = 0; m < curr.data.size(); m++) {
                    of<<curr.data.at(m)<<",";
                }
                of<<curr.t<<" to ";

                for(size_t m = 0; m < successor.data.size(); m++) {
                    of<<successor.data.at(m)<<",";
                }
                of<<successor.t<<" as 1\n";
            }
        }

        of << "splot '-' using 1:2:3 with points palette pointtype 7 ps 1; pause -1;\n";

        for(size_t i = 0; i<mst.size(); i++) {
            for(size_t m = 0; m < pc->at(0).data.size(); m++) {
                of << pc->at(i).data.at(m)<<" ";
            }
            of << pc->at(i).t <<"\n";
        }

        of.close();
        return true;
    }    

    /*
    bool DirectedMST::draw_speed(const char * fname) {
        std::ofstream of(fname);
        of << std::fixed;  // set float style
        if (!of.is_open()) {
            std::cerr << "[Error] could not save under '" << fname << std::endl;
            return false;
        }

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
        Point0 a;
        Point0 b;
        of <<"e\n";
        double scale = 10.0;
        for(size_t i = 0; i<pca_equations_n1.size(); i++) { 
            a = pca_equations_n1.at(i).at(0);
            b = pca_equations_n1.at(i).at(1);
            of<<a.x<<" "<<a.y<<" "<<a.t<<"\n";
            of<<a.x+scale*b.x<<" "<<a.y+scale*b.y<<" "<<a.t+scale*b.t<<"\n\n\n";
        }
        
        of <<"e\n";
        scale = 10.0;
        for(size_t i = 0; i<pca_equations_n2.size(); i++) { 
            a = pca_equations_n2.at(i).at(0);
            b = pca_equations_n2.at(i).at(1);
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
            Point0 curr; 
            Point0 successor;
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
            Point0 curr;
            Point0 successor;
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

*/