#ifndef POINTCLOUD_HPP
#define POINTCLOUD_HPP

#include <vector>
#include <ostream>
#include <numeric>

//#include <iostream>

class Point0 {
 public:
  std::vector<double> data;
  double t;

  Point0(){this->data.clear(); this->t=0.0;};
  Point0(size_t ndims){for(size_t i = 0; i < (ndims-1); i++) { this->data.push_back(0.0); }; this->t=0.0; };
  Point0(std::vector<double>& data, double t);

  Point0(const Point0& t) { 
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

  double euclidian_distance(const Point0& other);
  double no_root_euclidian_distance(const Point0& other);
  double spacial_euclidian_distance(const Point0& other);

  friend std::ostream& operator<<(std::ostream& os, const Point0& p);
  bool operator==(const Point0& p) const; 
  Point0& operator=(const Point0& other);
  // vector addition
  Point0 operator+(const Point0& p) const;
  // vector subtraction
  Point0 operator-(const Point0& p) const;
  // scalar product
  double operator*(const Point0& p) const;
  // scalar division
  Point0 operator/(double c) const;
  // scalar multiplication
  Point0 operator*(double c) const;
  double angle(const Point0 & p);
  double cosine_similarity(const Point0 & p);

};


// The Pointcloud is a vector of points
class PointCloud : public std::vector<Point0> {};


double orthogonal_lsq(PointCloud & pc, Point0* a, Point0* b);
double orthogonal_lsq_distance(Point0* a, Point0* b, Point0 * c);
float angle(Point0 * a, Point0 * b);


#endif
