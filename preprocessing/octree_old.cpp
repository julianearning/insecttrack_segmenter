#include <iostream>
#include "octree.hpp"


void choose_point() {
    
}


BoundingBox& BoundingBox::operator=(const BoundingBox& other) {
    upper_left.x = other.upper_left.x;
    upper_left.y = other.upper_left.y;
    upper_left.t = other.upper_left.t;
    dimension = other.dimension;
    return *this;
}


bool check_if_in_bb(Point & p, BoundingBox & bb) {
    return ((bb.upper_left.x <= p.x) && (bb.upper_left.y <= p.y) && (bb.upper_left.t <= p.t)) && ( ((bb.upper_left.x+bb.dimension) >= p.x) && ((bb.upper_left.y+bb.dimension)>=p.y) && ((bb.upper_left.t+bb.dimension)>=p.t));
}

OcTree::OcTree(PointCloud & pc, int max_depth, int curr_depth, BoundingBox & bb) {
    this->curr_depth = curr_depth;
    if(curr_depth>=max_depth) {
        return;
    }
    curr_depth++;
    this->m_objects = pc;
    Point p;
    p.x=bb.upper_left.x;
    p.y=bb.upper_left.y;
    p.t=bb.upper_left.t;
    this->m_region = BoundingBox(p, bb.dimension);
    this->max_depth = max_depth;
    this->_parent = NULL;
    int dimension_halved=bb.dimension/2;
    int start=bb.upper_left.x;
    int end=start+dimension_halved;
    std::cout<<m_region.upper_left<<", Dimension: "<<bb.dimension<<", "<<curr_depth<<"\n";
    //char combination=0;
    //char mask0=1;
    //char mask1=2;
    //char mask2=4;
//
    //for(int i = 0; i<8; i++) {   
    //    std::cout<<((mask2&combination)>>2)<<((mask1&combination)>>1)<<(mask0&combination)<<std::endl;
    //    combination++;
    //}


    
    // Box 0: start,start,start,dimension
    p.x=start;
    p.y=start;
    p.t=start;
    BoundingBox a = BoundingBox(p,dimension_halved);
    //std::cout<<"Box 0: "<<a.upper_left<<"\n";
    OcTree(pc,max_depth,curr_depth,a);
    // Box 1: end,start,start,dimension
    p.x=end;
    p.y=start;
    p.t=start;
    a = BoundingBox(p,dimension_halved);
    //std::cout<<"Box 1: "<<a.upper_left<<"\n";
    OcTree(pc,max_depth,curr_depth,a);
    // Box 2: start,end,start,dimension
    p.x=start;
    p.y=end;
    p.t=start;
    a = BoundingBox(p,dimension_halved);
    //std::cout<<"Box 2: "<<a.upper_left<<"\n";
    OcTree(pc,max_depth,curr_depth,a);
    // Box 3: start,start,end,dimension
    p.x=start;
    p.y=start;
    p.t=end;
    a = BoundingBox(p,dimension_halved);
    //std::cout<<"Box 3: "<<a.upper_left<<"\n";
    OcTree(pc,max_depth,curr_depth,a);
    // Box 4: end,end,start,dimension
    p.x=end;
    p.y=end;
    p.t=start;
    a = BoundingBox(p,dimension_halved);
    //std::cout<<"Box 4: "<<a.upper_left<<"\n";
    OcTree(pc,max_depth,curr_depth,a);
    // Box 5: end,start,end,dimension
    p.x=end;
    p.y=start;
    p.t=end;
    a = BoundingBox(p,dimension_halved);
    //std::cout<<"Box 5: "<<a.upper_left<<"\n";
    OcTree(pc,max_depth,curr_depth,a);
    // Box 6: start,end,end,dimension
    p.x=start;
    p.y=end;
    p.t=end;
    a = BoundingBox(p,dimension_halved);
    //std::cout<<"Box 6: "<<a.upper_left<<"\n";
    OcTree(pc,max_depth,curr_depth,a);
    // Box 7: end,end,end,dimension
    p.x=end;
    p.y=end;
    p.t=end;
    a = BoundingBox(p,dimension_halved);
    //std::cout<<"Box 7: "<<a.upper_left<<"\n";
    OcTree(pc,max_depth,curr_depth,a);
    //for (std::vector<Point>::iterator i = pc.begin(); i != pc.end(); ++i) {
    //    
    //}
    
}