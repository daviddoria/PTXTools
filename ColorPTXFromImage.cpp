// STL
#include <fstream>
#include <iostream>
#include <vector>

// Eigen
#include <Eigen/Geometry>
#include <Eigen/Dense> // for Vector

// VTK
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>

// ITK
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageRegionConstIterator.h"

// Custom
#include "PTXImage.h"

Eigen::MatrixXd ReadP(const std::string& filename);

int main(int argc, char *argv[])
{
  if(argc < 5)
    {
    std::cerr << "Required arguments: cameraMatrixFileName.txt ptxFileName.ptx imageFileName.png outputFileName.ptx" << std::endl;
    return EXIT_FAILURE;
    }
  // Parse arguments
  std::string cameraMatrixFileName = argv[1];
  std::string ptxFileName = argv[2];
  std::string imageFileName = argv[3];
  std::string outputFileName = argv[4];

  // Output arguments
  std::cout << "cameraMatrixFileName : " << cameraMatrixFileName << std::endl;
  std::cout << "ptxFileName : " << ptxFileName << std::endl;
  std::cout << "imageFileName : " << imageFileName << std::endl;
  std::cout << "outputFileName : " << outputFileName << std::endl;

  // Read the PTX file
  PTXImage ptxImage;
  ptxImage.ReadFile(ptxFileName);
  ptxImage.WritePTX("test.ptx");
  
  PTXImage::XYZImageType::Pointer xyzImage = ptxImage.GetXYZImage();
  
  PTXImage::RGBImageType::Pointer colorImage = PTXImage::RGBImageType::New();
  ptxImage.CreateRGBImage(colorImage);
  std::cout << "colorImage: " << colorImage->GetLargestPossibleRegion() << std::endl;

  // Read the camera matrix relating the input image to the PTX/scan/LiDAR file
  Eigen::MatrixXd P = ReadP(cameraMatrixFileName);
  std::cout << "P: " << P << std::endl;

  // Read the input image
  typedef itk::ImageFileReader<PTXImage::RGBImageType> ReaderType;
  ReaderType::Pointer imageReader = ReaderType::New();
  imageReader->SetFileName(imageFileName);
  imageReader->Update();

  itk::ImageRegionConstIterator<PTXImage::XYZImageType> xyzImageIterator(xyzImage, xyzImage->GetLargestPossibleRegion());
  //ptxImage.WritePTX("test2.ptx");
  while(!xyzImageIterator.IsAtEnd())
  {
    
    // Get the value of the current pixel
    PTXImage::XYZImageType::PixelType xyz = xyzImageIterator.Get();

    Eigen::Vector4d X;
    X(0) = xyz[0];
    X(1) = xyz[1];
    X(2) = xyz[2];
    X(3) = 1;

    // Perform the projection
    Eigen::Vector3d projectedHomog = P * X;

    // Get the projected point in non-homogeneous coordinates
    Eigen::Vector2d projected = projectedHomog.hnormalized();

    //std::cout << "Projected " << p[0] << " " << p[1] << " " << p[2] << " to : " << projected << std::endl;

    // Get the pixel that the point was projected to
    itk::Index<2> projectedPixel;
    projectedPixel[0] = round(projected(0));
    projectedPixel[1] = round(projected(1));

    //std::cout << "pixel: " << projectedPixel << std::endl;

    PTXImage::RGBImageType::PixelType color;
    
    if(!imageReader->GetOutput()->GetLargestPossibleRegion().IsInside(projectedPixel))
      {
      std::cout << "Point does not project to image!" << std::endl;
      color.SetRed(0);
      color.SetGreen(255);
      color.SetBlue(0);
      }
    else
      {
      color = imageReader->GetOutput()->GetPixel(projectedPixel);
      }
    
    //std::cout << "color: " << color << std::endl;

    colorImage->SetPixel(xyzImageIterator.GetIndex(), color);

    ++xyzImageIterator;
  }

  ptxImage.WritePTX("test3.ptx");

  ptxImage.ReplaceRGB(colorImage);
  ptxImage.WritePTX("test4.ptx");
  
  ptxImage.WritePointCloud("outputPointCloud");
  ptxImage.WritePTX("output.ptx");

  return EXIT_SUCCESS;
}

Eigen::MatrixXd ReadP(const std::string& filename)
{
  Eigen::MatrixXd P(3,4);

  std::string line;
  std::ifstream fin(filename.c_str());

  if(fin == NULL)
  {
    std::cout << "Cannot open file." << std::endl;
  }

  unsigned int lineCounter = 0;
  while(getline(fin, line))
    {
    std::stringstream ss;
    ss << line;
    ss >> P(lineCounter,0);
    ss >> P(lineCounter,1);
    ss >> P(lineCounter,2);
    ss >> P(lineCounter,3);
    lineCounter++;
    }
  return P;

}
