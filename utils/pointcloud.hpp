#ifndef POINTCLOUD_HPP
#define POINTCLOUD_HPP

#include <vector>
#include <ostream>
#include <numeric>
#include "../kdtree/kdtree.hpp"

//#include <iostream>

class Point {
 public:
  std::vector<double> data;
  double t;

  Point(){this->data.clear(); this->t=0.0;};
  Point(size_t ndims){for(size_t i = 0; i < (ndims-1); i++) { this->data.push_back(0.0); }; this->t=0.0; };
  Point(std::vector<double>& data, double t);

  Point(const Point& t) { 
    //std::cout<<"Copy constructor ;)\n";
    data.clear();
    for(size_t i = 0; i<t.data.size(); i++) {
      data.push_back(t.data.at(i));
    }
    this->t=t.t;
  }

  // representation of 3D point as std::vector
  std::vector<double> as_vector() const;
  // Euclidean norm
  double norm() const;
  // squared norm
  double unsquared_norm() const;

  double euclidian_distance(const Point& other);
  double no_root_euclidian_distance(const Point& other);
  double spacial_euclidian_distance(const Point& other);

  friend std::ostream& operator<<(std::ostream& os, const Point& p);
  bool operator==(const Point& p) const; 
  Point& operator=(const Point& other);
  // vector addition
  Point operator+(const Point& p) const;
  // vector subtraction
  Point operator-(const Point& p) const;
  // scalar product
  double operator*(const Point& p) const;
  // scalar division
  Point operator/(double c) const;
  // scalar multiplication
  Point operator*(double c) const;
  double angle(const Point & p);
  double cosine_similarity(const Point & p);

};


// The Pointcloud is a vector of points
class PointCloud : public std::vector<Point> {};


double orthogonal_lsq(PointCloud & pc, Point* a, Point* b);
double orthogonal_lsq_distance(Point* a, Point* b, Point * c);
float angle(Point * a, Point * b);


#endif
