#include <algorithm>
#include <cmath>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <cmath>
//#include <iostream>

#include "pointcloud.hpp"


Point0::Point0(std::vector<double>& data, double t) {
  this->data.clear();
  for(size_t i = 0; i<data.size(); i++) {
    this->data.push_back(data.at(i));
  }
  this->t = t;
}

// representation of 3D point as std::vector.
std::vector<double> Point0::as_vector() const {
  std::vector<double> point(this->data.size()+1);
  for(size_t i = 0; i<data.size(); i++) {
    point.push_back(this->data.at(i));
  }
  point.push_back(this->t);
  return point;
}

// Euclidean norm of the point
//double Point0::norm() const { return sqrt((x * x) + (y * y) + (t * t)); }
double Point0::norm() const { 
  double res=0.0;
  for(size_t i = 0; i<data.size(); i++) {
    res+=(data.at(i)*data.at(i));
  }
  return sqrt(res+(t*t)); 
}

double Point0::euclidian_distance(const Point0& other) {
  if(other.data.size()!=this->data.size()) return -1; 
  double res=0.0;
  double curr=0.0;
  for(size_t i=0; i<this->data.size(); i++) {
    curr=other.data.at(i)-this->data.at(i);
    res+=(curr*curr);
  }
  curr = other.t - this->t; 
  return std::sqrt( res + (curr * curr));
}

double Point0::spacial_euclidian_distance(const Point0& other) {
  if(other.data.size()!=this->data.size()) return -1; 
  double res=0.0;
  double curr=0.0;
  for(size_t i=0; i<this->data.size(); i++) {
    curr=other.data.at(i)-this->data.at(i);
    res+=(curr*curr);
  }
  return std::sqrt( res );
}

double Point0::no_root_euclidian_distance(const Point0& other) {
  if(other.data.size()!=this->data.size()) return -1; 
  double res=0.0;
  double curr=0.0;
  for(size_t i=0; i<this->data.size(); i++) {
    curr=other.data.at(i)-this->data.at(i);
    res+=(curr*curr);
  }
  curr = other.t - this->t; 
  return (res + (curr * curr));
}

double Point0::cosine_similarity(const Point0& p) {

  double res=0.0;
  for(size_t i = 0; i<this->data.size(); i++) {
    res+=(this->data.at(i)*p.data.at(i));
  }
  res+=(this->t*p.t);

  return res / (this->norm()*p.norm());
  //return (this->x*p.x+this->y*p.y+this->t*p.t) / (this->norm()*p.norm());
}

double Point0::angle(const Point0& p) {
  double res=0.0;
  for(size_t i = 0; i<this->data.size(); i++) {
    res+=(this->data.at(i)*p.data.at(i));
  }
  res+=(this->t*p.t);
  return std::acos(res/(this->norm()*p.norm()));
  //return std::acos((float)(this->x*p.x+this->y*p.y+this->t*p.t)/(std::sqrt(this->x*this->x+this->y*this->y+this->t*this->t)*std::sqrt(p.x*p.x+p.y*p.y+p.t*p.t)));
}


// squared euclidean norm of the point
//double Point0::unsquared_norm() const { return (x * x) + (y * y) + (t * t); }
double Point0::unsquared_norm() const { 
  double res=0.0;
  for(size_t i = 0; i<this->data.size(); i++) {
    res+=(this->data.at(i)*this->data.at(i));
  }
  res+=(this->t*this->t);
  return res;
}

Point0 &Point0::operator=(const Point0 &other) {
  //std::cout<<"Assignment ;)\n";
  this->data.clear();
  for(size_t i = 0; i<other.data.size(); i++) {
    this->data.push_back(other.data.at(i));
  }
  this->t = other.t;
  return *this;
}

bool Point0::operator==(const Point0 &p) const {
  //return (x == p.x && y == p.y && t == p.t);   // Rundungsfehler   
  //return ((x == p.x) && (y == p.y) && (std::abs(t-p.t) < 0.1));   // approx. version 
  for(size_t i = 0; i<this->data.size(); i++) {
    if(this->data.at(i)!=p.data.at(i)) return false;
  }
  return (std::abs(t-p.t) < 0.1);
}

// formatted output of the pointer
std::ostream &operator<<(std::ostream &strm, const Point0 &p) {
  for(size_t i = 0; i<p.data.size(); i++) {
    strm << p.data.at(i) << " ";
  }
  strm << p.t;
  return strm; 
  //return strm << p.x << " " << p.y << " " << p.t;
}

// vector addition
Point0 Point0::operator+(const Point0 &p) const {
  //Point0 ret;
  //ret.x = this->x + p.x;
  //ret.y = this->y + p.y;
  //ret.t = this->t + p.t;
  //return ret;
  Point0 ret;
  for(size_t i = 0; i<this->data.size(); i++) {
    ret.data.push_back(this->data.at(i)+p.data.at(i));
  }
  ret.t = this->t+p.t;
  return ret;
}

// vector subtraction
Point0 Point0::operator-(const Point0 &p) const {
  //Point0 ret;
  //ret.x = this->x - p.x;
  //ret.y = this->y - p.y;
  //ret.t = this->t - p.t;
  //return ret;
  Point0 ret;
  for(size_t i = 0; i<this->data.size(); i++) {
    ret.data.push_back(this->data.at(i)-p.data.at(i));
  }
  ret.t = this->t-p.t;
  return ret;
}

// scalar product (dot product)
double Point0::operator*(const Point0 &p) const {
  double res=0.0;
  for(size_t i = 0; i<this->data.size(); i++) {
    res+=(this->data.at(i)*p.data.at(i));
  }
  res+=(this->t*p.t);
  return res;
  //return this->x * p.x + this->y * p.y + this->t * p.t;
}

// scalar division
Point0 Point0::operator/(double c) const { 
  std::vector<double> data_new;
  for(size_t i = 0; i<this->data.size(); i++) {
    data_new.push_back(this->data.at(i)/c);
  }
  return Point0(data_new, this->t/c);
  //return Point0(x / c, y / c, t / c); 
}

// scalar multiplication
Point0 Point0::operator*(double c) const { 
  std::vector<double> data_new;
  for(size_t i = 0; i<this->data.size(); i++) {
    data_new.push_back(this->data.at(i)*c);
  }
  return Point0(data_new, this->t*c);
}

// scalar multiplication
/* Point0 Point0::operator*(double c) {
  Point0 v(c * this->x, c * this->y, c * this->t);
  return v;
} */

/* Point0 operator* (double c, const Point0& p) {
  Point0 v(c * p.x, c * p.y, c * p.t);
  return v; 
} */

