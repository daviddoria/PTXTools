#include "CameraCalibration.h"

#include <iostream>
#include <fstream>

Point2DVector LoadPoints2D(const std::string& filename)
{
  std::cout << "LoadPoint2D " << filename << std::endl;
  std::string line;
  std::ifstream fin(filename.c_str());
  Point2DVector points;
  if(fin == NULL)
    {
    std::cout << "Cannot open file." << std::endl;
    }

  while(getline(fin, line))
    {
    std::stringstream ss;
    ss << line;
    double p[3];
    ss >> p[0] >> p[1];
    points.push_back(Eigen::Vector2d (p[0], p[1]));
    }
  return points;
}

Point3DVector LoadPoints3D(const std::string& filename)
{
  std::cout << "LoadPoint3D " << filename << std::endl;
  std::string line;
  std::ifstream fin(filename.c_str());
  Point3DVector points;
  if(fin == NULL)
    {
    std::cout << "Cannot open file." << std::endl;
    }

  while(getline(fin, line))
    {
    std::stringstream ss;
    ss << line;
    double p[3];
    ss >> p[0] >> p[1] >> p[2];
    points.push_back(Eigen::Vector3d (p[0], p[1], p[2]));
    }
  return points;
}

Eigen::MatrixXd ComputeP_NormalizedDLT(const Point2DVector& points2D, const Point3DVector& points3D)
{
  unsigned int numberOfPoints = points2D.size();
  if(points3D.size() != numberOfPoints)
    {
    std::cerr << "The number of 2D points (" << points2D.size() << ") must match the number of 3D points (" << points3D.size() << ")!" << std::endl;
    exit(-1);
    }

  Eigen::MatrixXd similarityTransform2D = ComputeNormalizationTransform<Eigen::Vector2d>(points2D);
  Eigen::MatrixXd similarityTransform3D = ComputeNormalizationTransform<Eigen::Vector3d>(points3D);

//   std::cout << "Computed similarity transforms:" << std::endl;
//   std::cout << "similarityTransform2D: " << similarityTransform2D << std::endl;
//   std::cout << "similarityTransform3D: " << similarityTransform3D << std::endl;
//   
  Point2DVector transformed2DPoints(numberOfPoints);
  Point3DVector transformed3DPoints(numberOfPoints);

  for(unsigned int i = 0; i < numberOfPoints; ++i)
    {
    Eigen::VectorXd point2Dhomogeneous = points2D[i].homogeneous();
    Eigen::VectorXd point2Dtransformed = similarityTransform2D * point2Dhomogeneous;
    transformed2DPoints[i] = point2Dtransformed.hnormalized();

    Eigen::VectorXd point3Dhomogeneous = points3D[i].homogeneous();
    Eigen::VectorXd point3Dtransformed = similarityTransform3D * point3Dhomogeneous;
    transformed3DPoints[i] = point3Dtransformed.hnormalized();
  
    //transformed2DPoints[i] = (similarityTransform2D * points2D[i].homogeneous()).hnormalized();
    //transformed3DPoints[i] = (similarityTransform3D * points3D[i].homogeneous()).hnormalized();
    }

  std::cout << "Transformed points." << std::endl;
  
  // Compute the Camera Projection Matrix

  Eigen::MatrixXd A(2*numberOfPoints,12);
  for(unsigned int i = 0; i < numberOfPoints; ++i)
    {
    // First row/equation from the ith correspondence
    unsigned int row = 2*i;
    A(row, 0) = 0;
    A(row, 1) = 0;
    A(row, 2) = 0;
    A(row, 3) = 0;
    A(row, 4) = transformed3DPoints[i](0);
    A(row, 5) = transformed3DPoints[i](1);
    A(row, 6) = transformed3DPoints[i](2);
    A(row, 7) = 1;
    A(row, 8) = -transformed2DPoints[i](1) * transformed3DPoints[i](0);
    A(row, 9) = -transformed2DPoints[i](1) * transformed3DPoints[i](1);
    A(row, 10) = -transformed2DPoints[i](1) * transformed3DPoints[i](2);
    A(row, 11) = -transformed2DPoints[i](1);

    // Second row/equation from the ith correspondence
    row = 2*i+1;
    A(row, 0) = transformed3DPoints[i](0);
    A(row, 1) = transformed3DPoints[i](1);
    A(row, 2) = transformed3DPoints[i](2);
    A(row, 3) = 1;
    A(row, 4) = 0;
    A(row, 5) = 0;
    A(row, 6) = 0;
    A(row, 7) = 0;
    A(row, 8) = -transformed2DPoints[i](0) * transformed3DPoints[i](0);
    A(row, 9) = -transformed2DPoints[i](0) * transformed3DPoints[i](1);
    A(row, 10) = -transformed2DPoints[i](0) * transformed3DPoints[i](2);
    A(row, 11) = -transformed2DPoints[i](0);
    }

  std::cout << "A: " << A << std::endl;
  
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);

  Eigen::MatrixXd V = svd.matrixV();
  Eigen::MatrixXd lastColumnOfV = V.col(11);

  Eigen::MatrixXd P = Reshape(lastColumnOfV, 3, 4);
    
  // Denormalization
  P = similarityTransform2D.inverse()*P*similarityTransform3D; // 3x3 * 3x4 * 4x4 = 4x4

  return P;
}

Eigen::MatrixXd ComputeP_Nonlinear(const Point2DVector& points2d, const Point3DVector& points3d)
{
  Eigen::MatrixXd linearP = ComputeP_NormalizedDLT(points2d, points3d);
  return linearP;
}


Eigen::MatrixXd Reshape(const Eigen::VectorXd& vec, const unsigned int rows, const unsigned int cols)
{
  if(static_cast<unsigned int>(vec.rows()) != rows*cols)
    {
    std::cerr << "Cannot reshape a vector with " << vec.rows() << " to a " << rows << " x " << cols << " matrix!" << std::endl;
    exit(-1);
    }

  Eigen::MatrixXd P(rows,cols);
  for(unsigned int row = 0; row < rows; ++row)
    {
    for(unsigned int col = 0; col < cols; ++col)
      {
      P(row, col) = vec(row*cols + col);
      }
    }

  return P;
}


struct LMFunctor
{
  int operator()(const Eigen::VectorXf &x, Eigen::VectorXf &fvec) const
  {
    // Implement y = (x-5)^2 (remember, operator() should return the value BEFORE it is squared.
    fvec(0) = x(0) - 5.0;
    return 0;
  }

  int df(const Eigen::VectorXf &x, Eigen::MatrixXf &fjac) const
  {
    Eigen::VectorXf epsilon(1);
    epsilon(0) = 1e-5;

    Eigen::VectorXf fvec1(1);
    operator()(x + epsilon, fvec1);
    Eigen::VectorXf fvec2(1);
    operator()(x - epsilon, fvec2);
    fjac = (fvec1 - fvec2)/2.0f;
    return 0;
  }

  int inputs() const { return 1; }// inputs is the dimension of x.
  int values() const { return 1; } // "values" is the number of f_i and

  Point2DVector& points2d;
  Point3DVector& points3d;
};
/*
float NonLinearProjectionError(Eigen::Vector2d& parameters)
{
  
//   Parameters must be
//   parameters = [r(0) r(1) r(2) t(0) t(1) t(2) K(0,0), K(1,1), K(0,2), K(1,2), d]
//                [ 0    1    2    3    4    5    6        7       8       9     10
//   where r() is the vector of rodrigues angles, t is the translation vector, K is the intrinsic camera matrix, and d is the ?
  
  
  
  Eigen::MatrixXd R = rodrigues(k(1:3)');
  Eigen::Vector3d t;
  t(0) = parameters(3);
  t(1) = parameters(4);
  t(2) = parameters(5);
  
  float d = parameters(10);
  
  Eigen::MatrixXd K(3,3);
  K(0,0) = parameters(6);
  K(0,1) = 0;
  K(0,2) = parameters(8);
  K(1,0) = 0;
  K(1,1) = parameters(7);
  K(1,2) = parameters(9);
  K(0,0) = 0;
  K(2,1) = 0;
  K(2,2) = 1;

  for(unsigned int i = 0; i < points2d.size(); ++i)
    {
    Xw=X(1:3,i);
    Xc=R*Xw+t;

    xid=K*Xc;

    xid=xid/xid(3);

    xd=(xid(1)-K(1,3))/K(1,1);
    yd=(xid(2)-K(2,3))/K(2,2);



    r=xd^2+yd^2;

    xu=xd*(1+d*r);
    yu=yd*(1+d*r);

    xi(1)=xu*K(1,1)+K(1,3);
    xi(2)=yu*K(2,2)+K(2,3);



    err(i)=norm(xi'-x(1:2,i));
    }
}
*/
