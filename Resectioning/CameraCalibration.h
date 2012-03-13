#ifndef CameraCalibration_H
#define CameraCalibration_H

#include <vector>

#include <Eigen/Dense>

namespace CameraCalibration
{

Eigen::VectorXd GetCameraCenter(const Eigen::MatrixXd& P);

typedef std::vector<Eigen::Vector2d,Eigen::aligned_allocator<Eigen::Vector2d> > Point2DVector;
typedef std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> > Point3DVector;

/** This function reads 2D points from a file. */
Point2DVector LoadPoints2D(const std::string& filename);

/** This function reads 3D points from a file. */
Point3DVector LoadPoints3D(const std::string& filename);

/** This function computes the Normalized DLT solution of P from corresponding 2D and 3D points. */
Eigen::MatrixXd ComputeP_NormalizedDLT(const Point2DVector&, const Point3DVector&);

/** This function performs a nonlinear improvement to the DLT estimation.
 * TODO: Make sure this works (I don't think it does currently). */
Eigen::MatrixXd ComputeP_Nonlinear(const Point2DVector&, const Point3DVector&);

template<typename T>
T Centroid(const typename std::vector<T,typename Eigen::aligned_allocator<T> >& points);

template<typename T>
Eigen::MatrixXd ComputeNormalizationTransform(const typename
                std::vector<T,typename Eigen::aligned_allocator<T> >& points);

Eigen::MatrixXd HomogeneousMultiply(const Point3DVector& points);

Eigen::MatrixXd Reshape(const Eigen::VectorXd& vec, const unsigned int rows, const unsigned int cols);

}

#include "CameraCalibration.hpp"

#endif
