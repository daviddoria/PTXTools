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
  PTXImage ResectionSmart(const Eigen::MatrixXd& P, const PTXImage& ptxImage,
                          const PTXImage::RGBImageType* const colorImage);

  Eigen::MatrixXd ReadP(const std::string& filename);
  
  PTXImage ResectionNaive(const Eigen::MatrixXd& P, const PTXImage& ptxImage,
                          const PTXImage::RGBImageType* const colorImage);
}

#endif
