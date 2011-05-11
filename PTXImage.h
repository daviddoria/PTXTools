#ifndef PTXImage_h
#define PTXImage_h

// ITK
#include "itkImage.h"
#include "itkPoint.h"
#include "itkImageRegionIterator.h"

// VTK
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>

#include "PTXPixel.h"

/** \class PTXImage
 *  \brief This class handles common operations (reading, writing, appending, etc) operations on Leica PTX files.
 */
class PTXImage
{
public:

  // The image which contains complete information about each point
  typedef itk::Image<PTXPixel, 2> FullImageType;

  // A normal RGB image
  typedef itk::Image<itk::CovariantVector<unsigned char, 3>, 2> RGBImageType;

  // A binary image
  typedef itk::Image<unsigned char, 2> MaskImageType;

  // A scalar float image
  typedef itk::Image<float, 2> FloatImageType;
  
  // An image to hold X,Y,and Z coordinate channels.
  typedef itk::Image<itk::CovariantVector<float, 3>, 2> XYZImageType;

  // A 4 channel image (R, G, B, Depth)
  typedef itk::Image<itk::CovariantVector<float, 4>, 2> RGBDImageType;

  // A 5 channel image (R, G, B, Depth, Intensity)
  typedef itk::Image<itk::CovariantVector<float, 5>, 2> RGBDIImageType;

  // Constructor
  PTXImage();

  // Add/append a ptx file to the right of this one.
  void AppendPTXRight(PTXImage);
  
  // Coordinate images
  FloatImageType::Pointer GetCoordinateImage(unsigned int dimension);
  
  void WriteXYZ(const std::string filePrefix);
  XYZImageType::Pointer GetXYZImage();
  
  void WriteX(const std::string filePrefix);
  FloatImageType::Pointer GetXImage();
  
  void WriteY(const std::string filePrefix);
  FloatImageType::Pointer GetYImage();
  
  void WriteZ(const std::string filePrefix);
  FloatImageType::Pointer GetZImage();
  
  void WriteFloatImage(FloatImageType::Pointer image, const std::string filename);
  
  // Derivatives
  FloatImageType::Pointer GetLaplacian(unsigned int dimension);
  
  void WriteXYZLaplacian(const std::string filePrefix);
  XYZImageType::Pointer GetXYZLaplacian();
  
  void WriteXLaplacian(const std::string filePrefix);
  FloatImageType::Pointer GetXLaplacian();
  
  void WriteYLaplacian(const std::string filePrefix);
  FloatImageType::Pointer GetYLaplacian();
  
  void WriteZLaplacian(const std::string filePrefix);
  FloatImageType::Pointer GetZLaplacian();

  // Set the size of the PTXImage (this is done automatically if you use ReadFile)
  void SetSize(itk::ImageRegion<2>);
  
  // Downsample the ptx image by a factor of 'factor'
  PTXImage Downsample(const unsigned int factor);

  // Write a FullImageType to a ptx file
  void WritePTX(const std::string filename);

  // Write the projection plane and principal axis to a vtp file
  void WriteProjectionPlane(const std::string filename);

  void Crop(const itk::ImageRegion<2> region);

  // Get center direction ("principal axis")
  typedef itk::CovariantVector<float, 3> VectorType;
  VectorType GetPrincipalAxis();

  // Compute the Laplacian of the depth image where pixels are weighted by their distance to the center of the kernel.
  void ComputeWeightedDepthLaplacian(const std::string filename);

  // Create a 2D image of the points in the grid in which the points were acquired
  void CreateRGBImage(RGBImageType::Pointer image);
  void WriteRGBImage(std::string filename);

  // Create a colored point cloud
  void CreatePointCloud(vtkSmartPointer<vtkPolyData> pointCloud);
  void WritePointCloud(std::string outputFilePrefix);

  // Create an image of the intensities of the points in the grid in which they were acquired
  void CreateIntensityImage(FloatImageType::Pointer image);
  void WriteIntensityImage(std::string filePrefix);

  // Create a 2D, 1 channel image of the depths of the points in the grid in which they were acquired
  void CreateDepthImage(FloatImageType::Pointer image);
  void WriteDepthImage(std::string filePrefix);

  void WriteDepthLaplacian(const std::string filePrefix);
  FloatImageType::Pointer GetDepthLaplacian();
  
  void CreateRGBDImage(RGBDImageType::Pointer image);
  void WriteRGBDImage(std::string filePrefix);

  void CreateRGBDIImage(RGBDIImageType::Pointer image);
  void WriteRGBDIImage(std::string filePrefix);
  
  void WriteEverything(std::string filePrefix);

  // Actually read the PTX file
  void ReadFile(std::string filename);

  // This function allows the depth map to be modified externally and the new map applied to the grid
  void ReplaceDepth(itk::Image<float, 2>::Pointer depthImage);

  // This function allows the color and depth to be modified externally and the new map applied to the grid
  void ReplaceRGBD(itk::Image<itk::CovariantVector<float, 4>, 2>::Pointer rgbd);
  
  // This function allows the color to be modified
  void ReplaceRGB(itk::Image<itk::CovariantVector<float, 3>, 2>::Pointer rgb);
  
  // This function allows the color to be modified
  void ReplaceXYZ(itk::Image<itk::CovariantVector<float, 3>, 2>::Pointer xyz);

  // Blank the PTX image in areas where mask is non-zero
  void ApplyMask(itk::Image<unsigned char, 2>::Pointer mask);

  // Create a mask image where invalid pixels are non-zero
  void WriteInvalidMask(std::string& filename);

  // Create a mask image of points below a certain depthThreshold
  void WriteDepthThresholdMask(std::string& filename, float depthThreshold);

  // Count invalid points
  void CountInvalidPoints();

  // The main storage image.
  FullImageType::Pointer FullImage;

  // Access the main data.
  FullImageType::Pointer GetFullImage();

  // Get the size of the main data.
  itk::Size<2> GetSize();

  // Get the height of the main data.
  unsigned int GetHeight();

  // Get the width of the main data.
  unsigned int GetWidth();
  
  // Set a specific pixel to a specified value.
  void SetPixel(itk::Index<2>, PTXPixel&);

  // Find the nearest pixel which is not marked as invalid.
  itk::Index<2> FindNearestValidPixel(itk::Index<2> pixel, itk::Offset<2> offset);

  // Get the theta (side to side) angle of a pixel
  float ApproximateTheta(itk::Index<2> pixel);

  // Get the phi (up and down) angle of a pixel
  float ApproximatePhi(itk::Index<2> pixel);
  
  itk::Point<float, 3> ApproximateOldPoint(itk::Index<2> pixel);

  // Get the phi (up and down) angle of a pixel
  float GetPhi(itk::Index<2> index);

  // Get the theta (side to side) angle of a pixel
  float GetTheta(itk::Index<2> index);

  void ComputeAverageDeltaPhi();
  void ComputeAverageDeltaTheta();
  float AverageDeltaTheta;
  float AverageDeltaPhi;

  float DistanceBetweenPoints(PTXPixel a, PTXPixel b);

  PTXImage OrthogonalProjection(VectorType axis);

  itk::Index<2> FindValidTopCenterPixel();
  itk::Index<2> FindValidCenterPixel();
};


#endif