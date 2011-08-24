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
#include "Helpers.h"
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
  //FilePrefix prefix("test");
  //ptxImage.WritePTX(prefix);
  
  PTXImage::XYZImageType::Pointer xyzImage = ptxImage.GetXYZImage();
  
  PTXImage::RGBImageType::Pointer colorImage = PTXImage::RGBImageType::New();
  colorImage->SetRegions(xyzImage->GetLargestPossibleRegion());
  colorImage->Allocate();
  
  PTXImage::RGBImageType::PixelType green;
  green.SetRed(0);
  green.SetGreen(255);
  green.SetBlue(0);
  
  colorImage->FillBuffer(green);
  //ptxImage.CreateRGBImage(colorImage);

  //std::cout << "colorImage: " << colorImage->GetLargestPossibleRegion() << std::endl;

  // Read the camera matrix relating the input image to the PTX/scan/LiDAR file
  Eigen::MatrixXd P = ReadP(cameraMatrixFileName);
  std::cout << "P: " << P << std::endl;

  // Read the input image
  typedef itk::ImageFileReader<PTXImage::RGBImageType> ReaderType;
  ReaderType::Pointer imageReader = ReaderType::New();
  imageReader->SetFileName(imageFileName);
  imageReader->Update();

  typedef itk::Image<std::vector<itk::Index<2> >, 2> PixelImageType;
  PixelImageType::Pointer projectedImage = PixelImageType::New();
  projectedImage->SetRegions(imageReader->GetOutput()->GetLargestPossibleRegion());
  projectedImage->Allocate();

  std::vector<itk::Index<2> > emptyVector;
  //projectedImage->FillBuffer(emptyVector);
  Helpers::SetAllPixels<PixelImageType>(projectedImage, emptyVector);
  
  itk::ImageRegionConstIterator<PTXImage::XYZImageType> xyzImageIterator(xyzImage, xyzImage->GetLargestPossibleRegion());
  
  unsigned int badPoints = 0;

  // Iterate over the scan points image and track where each projects in the 'projectedImage'
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

    if(!imageReader->GetOutput()->GetLargestPossibleRegion().IsInside(projectedPixel))
      {
      // Do nothing
      }
    else
      {
      std::vector<itk::Index<2> > projectedSoFar = projectedImage->GetPixel(projectedPixel);
      projectedSoFar.push_back(xyzImageIterator.GetIndex());
      projectedImage->SetPixel(projectedPixel, projectedSoFar);
      //std::cout << "pixel: " << projectedPixel << std::endl;
      }

    ++xyzImageIterator;
    }

  std::cout << "There were " << badPoints << " that did not project to inside of the image!" << std::endl;

  itk::ImageRegionConstIterator<PixelImageType> projectedImageIterator(projectedImage, projectedImage->GetLargestPossibleRegion());

  while(!projectedImageIterator.IsAtEnd())
    {
    if(projectedImageIterator.Get().size() == 0)
      {
      //std::cout << "Point does not project to image!" << std::endl;
      badPoints++;
      }
    else
      {
      /*
      // Find and color the closest point
      itk::Index<2> closestPixel;
      float minDepth = std::numeric_limits<float>::max();
      //std::cout << "There were " << projectedImageIterator.Get().size() << " points that projected to this pixel." << std::endl;
      for(unsigned int i = 0; i < projectedImageIterator.Get().size(); ++i)
        {
        itk::Index<2> currentPixel = projectedImageIterator.Get()[i];
        //std::cout << "Current pixel: " << currentPixel << std::endl;

        // Determine the minimum depth
        PTXPixel ptxPixel = ptxImage.GetPTXPixel(currentPixel);
        if(ptxPixel.GetDepth() < minDepth)
          {
          minDepth = ptxPixel.GetDepth();
          closestPixel = projectedImageIterator.Get()[i];
          }
        }
      PTXImage::RGBImageType::PixelType color;
      color = imageReader->GetOutput()->GetPixel(projectedImageIterator.GetIndex());
      colorImage->SetPixel(closestPixel, color);
      */


      // Find the closest point
      float minDepth = std::numeric_limits<float>::max();
      //std::cout << "There were " << projectedImageIterator.Get().size() << " points that projected to this pixel." << std::endl;
      for(unsigned int i = 0; i < projectedImageIterator.Get().size(); ++i)
        {
        itk::Index<2> currentPixel = projectedImageIterator.Get()[i];
        //std::cout << "Current pixel: " << currentPixel << std::endl;

        // Determine the minimum depth
        PTXPixel ptxPixel = ptxImage.GetPTXPixel(currentPixel);
        if(ptxPixel.GetDepth() < minDepth)
          {
          minDepth = ptxPixel.GetDepth();
          }
        }
        
      // Color all points within a tolerance of the minimum depth
      PTXImage::RGBImageType::PixelType color;
      color = imageReader->GetOutput()->GetPixel(projectedImageIterator.GetIndex());
      for(unsigned int i = 0; i < projectedImageIterator.Get().size(); ++i)
        {
        itk::Index<2> currentPixel = projectedImageIterator.Get()[i];
        PTXPixel ptxPixel = ptxImage.GetPTXPixel(currentPixel);

        if(fabs(ptxPixel.GetDepth() - minDepth) < .1)
          {
          colorImage->SetPixel(currentPixel, color);
          }
        }

      } // end else over > 0 projections

    ++projectedImageIterator;
    } // end while over whole image
  
  ptxImage.ReplaceRGB(colorImage);
  
  //FilePrefix vtpPrefix("outputPointCloud");
  //ptxImage.WritePointCloud(vtpPrefix);
  
  FilePrefix ptxPrefix(outputFileName);
  ptxImage.WritePTX(ptxPrefix);

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
