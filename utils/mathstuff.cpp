#include "mathstuff.hpp"
#include <Eigen/Dense>
#include <iostream>

//
// orthogonal least squares fit with libeigen
//   pc: points
//   a, b: line representation as a + t*b
//   RC: largest eigenvalue
//
double MathStuff::orthogonal_lsq(PointCloud & pc, Point* a, Point* b, size_t ndims) {
    double rc = 0.0;

    if (pc.size() == 0)
      return rc;

    // anchor point is mean value
    //a->x = a->y = a->t = 0.0;
    *a = Point(ndims); 
    for (size_t i=0; i < pc.size(); i++) {
      *a = *a + pc.at(i);
    }
    *a = *a / pc.size();

    // copy points to libeigen matrix
    Eigen::MatrixXf points = Eigen::MatrixXf::Constant(pc.size(), ndims, 0);
    for (unsigned int i = 0; i < points.rows(); i++) {
      for(size_t k = 0; k < (ndims-1); k++) {
        points(i,k) = pc.at(i).data.at(k);
      }
      points(i,(ndims-1)) = pc.at(i).t;
    }   
    // compute scatter matrix ...
    Eigen::MatrixXf centered = points.rowwise() - points.colwise().mean();
    Eigen::MatrixXf scatter = (centered.adjoint() * centered); 
    // ... and its eigenvalues and eigenvectors
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> eig(scatter);
    Eigen::MatrixXf eigvecs = eig.eigenvectors();   
    // we need eigenvector to largest eigenvalue
    // libeigen yields it as LAST column
    //b->x = eigvecs(0,2); b->y = eigvecs(1,2); b->t = eigvecs(2,2);
    for(size_t k = 0; k < (ndims-1); k++) {
      b->data.at(k) = eigvecs(k,2);
    }
    b->t = eigvecs(ndims-1,2);
    //points(i,(ndims-1)) = pc.at(i).t;
    rc = eig.eigenvalues()(2);  
    return (rc);
}

double MathStuff::orthogonal_lsq_distance(Point* a, Point* b, Point * c, size_t ndims) {
  double distance = 0;
  double lambda;

  // copy points to libeigen vektor
  //Eigen::Vector3f a_vec = Eigen::Vector3f(a->x,a->y,a->t);
  //Eigen::Vector3f b_vec = Eigen::Vector3f(b->x,b->y,b->t);
  //Eigen::Vector3f c_vec = Eigen::Vector3f(c->x,c->y,c->t);

  Eigen::Vector3f a_vec = Eigen::Vector3f(ndims);
  Eigen::Vector3f b_vec = Eigen::Vector3f(ndims);
  Eigen::Vector3f c_vec = Eigen::Vector3f(ndims);

  for(size_t k = 0; k<(ndims-1); k++) {
    a_vec(k) = a->data.at(k);
    b_vec(k) = b->data.at(k);
    c_vec(k) = c->data.at(k);
  }
  a_vec(ndims-1) = a->t;
  b_vec(ndims-1) = b->t;
  c_vec(ndims-1) = c->t;


  //lambda = b_vec.transpose() * (c_vec - a_vec);
//
  //distance = (c_vec - (a_vec + lambda * b_vec)).norm();

  distance = (b_vec.cross(a_vec-c_vec)).norm()/b_vec.norm();



  return distance; 
}

