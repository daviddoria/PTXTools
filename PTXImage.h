/*
Copyright (C) 2011 David Doria, daviddoria@gmail.com

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef PTXImage_h
#define PTXImage_h

// ITK
#include "itkImage.h"
#include "itkPoint.h"
#include "itkImageRegionIterator.h"
#include "itkRGBPixel.h"

// VTK
#include <vtkSmartPointer.h>
class vtkPolyData;

// Custom
#include "PTXPixel.h"
#include "FilePrefix.h"

/** \class PTXImage
 *  \brief This class handles common operations (reading, writing, appending, etc) operations on Leica PTX files.
 */
class PTXImage
{
public:

  // The image which contains complete information about each point
  typedef itk::Image<PTXPixel, 2> FullImageType;

  // A normal RGB image
  typedef itk::Image<itk::RGBPixel<unsigned char>, 2> RGBImageType;
  typedef itk::Image<itk::CovariantVector<unsigned char, 3>, 2> RGBVectorImageType;

  // A binary image
  typedef itk::Image<unsigned char, 2> MaskImageType;
  typedef itk::Image<unsigned char, 2> UnsignedCharImageType;

  // A scalar float image
  typedef itk::Image<float, 2> FloatImageType;
  typedef itk::Image<float, 2> DepthImageType;

  // An image to hold X,Y,and Z coordinate channels.
  typedef itk::Image<itk::CovariantVector<float, 3>, 2> XYZImageType;

  // A 4 channel image (R, G, B, Depth)
  typedef itk::Image<itk::CovariantVector<float, 4>, 2> RGBDImageType;

  // A 4 channel image (R, G, B, Depth, Valid)
  typedef itk::Image<itk::CovariantVector<float, 5>, 2> RGBDVImageType;

  // A 5 channel image (R, G, B, Depth, Intensity)
  typedef itk::Image<itk::CovariantVector<float, 5>, 2> RGBDIImageType;

  // Constructor
  PTXImage();

  // Add/append a ptx file to the right of this one.
  void AppendPTXRight(const PTXImage&);

  // Coordinate images
  FloatImageType::Pointer GetCoordinateImage(const unsigned int dimension) const;

  void WriteXYZ(const std::string& filePrefix) const;
  XYZImageType::Pointer GetXYZImage() const;

  void WriteX(const std::string& filePrefix) const;
  FloatImageType::Pointer GetXImage() const;

  void WriteY(const std::string& filePrefix) const;
  FloatImageType::Pointer GetYImage() const;

  void WriteZ(const std::string& filePrefix) const;
  FloatImageType::Pointer GetZImage() const;

  void WriteFloatImage(const FloatImageType::Pointer image, const std::string& filename) const;

  // Derivatives
  FloatImageType::Pointer GetLaplacian(const unsigned int dimension) const;

  void WriteXYZLaplacian(const std::string& filePrefix) const;
  XYZImageType::Pointer GetXYZLaplacian() const;

  void WriteXLaplacian(const std::string& filePrefix) const;
  FloatImageType::Pointer GetXLaplacian() const;

  void WriteYLaplacian(const std::string& filePrefix) const;
  FloatImageType::Pointer GetYLaplacian() const;

  void WriteZLaplacian(const std::string& filePrefix) const;
  FloatImageType::Pointer GetZLaplacian() const;

  // Set the size of the PTXImage (this is done automatically if you use ReadFile)
  void SetSize(const itk::ImageRegion<2>&);

  // Downsample the ptx image by a factor of 'factor'
  PTXImage Downsample(const unsigned int factor) const;

  // Write a FullImageType to a ptx file
  void WritePTX(const FilePrefix& filename) const;

  // Write the projection plane and principal axis to a vtp file
  void WriteProjectionPlane(const std::string& filename) const;

  void Crop(const itk::ImageRegion<2>& region);

  PTXPixel GetPTXPixel(const itk::Index<2>& pixel) const;

  // Get center direction ("principal axis")
  typedef itk::CovariantVector<float, 3> VectorType;
  VectorType GetPrincipalAxis() const;

  // Compute the Laplacian of the depth image where pixels are weighted by their distance to the center of the kernel.
  void ComputeWeightedDepthLaplacian(const std::string& filename) const;

  // Create a 2D image of the points in the grid in which the points were acquired
  void CreateRGBImage(RGBImageType::Pointer image) const;
  void WriteRGBImage(const FilePrefix& prefix) const;

  // Create a colored point cloud
  void CreatePointCloud(vtkSmartPointer<vtkPolyData> pointCloud) const;
  void WritePointCloud(const FilePrefix& prefix) const;
  void WritePointCloud(const std::string& fileName) const;

  // Create an image of the intensities of the points in the grid in which they were acquired
  void CreateIntensityImage(FloatImageType::Pointer image) const;
  void WriteIntensityImage(const FilePrefix& filePrefix) const;

  // Create a 2D, 1 channel image of the depths of the points in the grid in which they were acquired
  void CreateDepthImage(FloatImageType::Pointer image) const;
  void WriteDepthImage(const FilePrefix& filePrefix) const;

  void WriteDepthLaplacian(const FilePrefix& filePrefix) const;
  FloatImageType::Pointer GetDepthLaplacian() const;

  void CreateRGBDImage(RGBDImageType::Pointer image) const;
  void WriteRGBDImage(const FilePrefix& filePrefix) const;

  void CreateRGBDVImage(RGBDVImageType::Pointer image) const;
  void WriteRGBDVImage(const FilePrefix& filePrefix) const;

  void CreateRGBDIImage(RGBDIImageType::Pointer image) const;
  void WriteRGBDIImage(const FilePrefix& filePrefix) const;

  void WriteEverything(const FilePrefix& filePrefix) const;

  // Actually read the PTX file
  void ReadFile(const std::string& filename);

  // This function allows the validity image to be modified externally and the new image applied to the grid
  void ReplaceValidity(const MaskImageType::Pointer validityImage);

  // This function allows the depth map to be modified externally and the new map applied to the grid
  void ReplaceDepth(const FloatImageType::Pointer depthImage);

  // This function allows the color and depth to be modified externally and the new map applied to the grid
  void ReplaceRGBD(const RGBDImageType::Pointer rgbd);

  // This function allows the color to be modified
  void ReplaceRGB(const RGBVectorImageType::Pointer rgb);
  void ReplaceRGB(const RGBImageType::Pointer rgb);

  // This function allows the color to be modified
  void ReplaceXYZ(const XYZImageType::Pointer xyz);

  // Blank the PTX image in areas where mask is non-zero
  void ApplyMask(const MaskImageType::Pointer mask);

  // Create a mask image where invalid pixels are non-zero
  void WriteInvalidMask(const std::string& filename) const;
  void CreateValidityImage(MaskImageType::Pointer image) const;

  // Create a mask image of points below a certain depthThreshold
  void WriteDepthThresholdMask(const std::string& filename, const float depthThreshold) const;

  // Count invalid points
  unsigned int CountInvalidPoints() const;
  unsigned int CountValidPoints() const;

  // The main storage image.
  FullImageType::Pointer FullImage;

  // Access the main data.
  FullImageType::Pointer GetFullImage() const;

  // Get the size of the main data.
  itk::Size<2> GetSize() const;

  // Get the height of the main data.
  unsigned int GetHeight() const;

  // Get the width of the main data.
  unsigned int GetWidth() const;

  // Set a specific pixel to a specified value.
  void SetPixel(const itk::Index<2>&, const PTXPixel&);

  // Find the nearest pixel which is not marked as invalid.
  itk::Index<2> FindNearestValidPixel(const itk::Index<2>& pixel, const itk::Offset<2>& offset) const;

  // Get the theta (side to side) angle of a pixel
  float ApproximateTheta(const itk::Index<2>& pixel) const;

  // Get the phi (up and down) angle of a pixel
  float ApproximatePhi(const itk::Index<2>& pixel) const;

  itk::Point<float, 3> ApproximateOldPoint(const itk::Index<2>& pixel) const;

  // Get the phi (up and down) angle of a pixel
  float GetPhi(const itk::Index<2>& index) const;

  // Get the theta (side to side) angle of a pixel
  float GetTheta(const itk::Index<2>& index) const;

  void ComputeAverageDeltaPhi();
  void ComputeAverageDeltaTheta();

  float DistanceBetweenPoints(const PTXPixel& a, const PTXPixel& b) const;

  PTXImage OrthogonalProjection(const VectorType& axis) const;

  itk::Index<2> FindValidTopCenterPixel() const;
  itk::Index<2> FindValidCenterPixel() const;

  void SetDebug(bool);

private:
  float AverageDeltaTheta;
  float AverageDeltaPhi;

  bool Debug;
};

#endif
