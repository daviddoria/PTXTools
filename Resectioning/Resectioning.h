#ifndef Resectioning_H
#define Resectioning_H

// Eigen
#include <Eigen/Dense> // for MatrixXd

// Custom
#include "PTXImage.h"

// STL
#include <string>

namespace Resectioning
{
  void ResectionSmart(const Eigen::MatrixXd& P, PTXImage& ptxImage, PTXImage::RGBImageType* const colorImage);

  Eigen::MatrixXd ReadP(const std::string& filename);
  
  // void ResectionNaive(Eigen::MatrixXd P, PTXImage& ptxImage, PTXImage::RGBImageType* const colorImage);
}

#endif
