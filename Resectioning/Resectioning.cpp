#include "Resectioning.h"

// Custom
#include "Helpers.h"
#include "ResectioningHelpers.h"
#include "PTXImage.h"

// Eigen
#include <Eigen/Geometry>
#include <Eigen/Dense> // for Vector

// STL
#include <iostream>
#include <string>

/** Color the provided PTX after mapping the colors from 'colorImage' through P. */
void ResectionSmart(Eigen::MatrixXd P, PTXImage& ptxImage, PTXImage::RGBImageType* const colorImage)
{
  std::cout << "Input has " << ptxImage.CountValidPoints() << " valid points." << std::endl;
  std::cout << "P: " << P << std::endl;

  //FilePrefix prefix("test");
  //ptxImage.WritePTX(prefix);

  PTXImage::XYZImageType::Pointer xyzImage = ptxImage.GetXYZImage();

  // This is the output image. Start by making it entirely green, then we will fill in valid values.
  PTXImage::RGBImageType::Pointer resultImage = PTXImage::RGBImageType::New();
  resultImage->SetRegions(xyzImage->GetLargestPossibleRegion());
  resultImage->Allocate();

  PTXImage::RGBImageType::PixelType green;
  green.SetRed(0);
  green.SetGreen(255);
  green.SetBlue(0);

  colorImage->FillBuffer(green);

  // We also want to track which pixels were filled
  PTXImage::MaskImageType::Pointer validityMask = PTXImage::MaskImageType::New();
  validityMask->SetRegions(xyzImage->GetLargestPossibleRegion());
  validityMask->Allocate();
  validityMask->FillBuffer(0);

  //ptxImage.CreateRGBImage(resultImage);

  //std::cout << "resultImage: " << resultImage->GetLargestPossibleRegion() << std::endl;

  typedef itk::Image<std::vector<itk::Index<2> >, 2> PixelImageType;
  PixelImageType::Pointer projectedImage = PixelImageType::New();
  projectedImage->SetRegions(ptxImage.GetFullRegion());
  projectedImage->Allocate();

  std::vector<itk::Index<2> > emptyVector;
  //projectedImage->FillBuffer(emptyVector);
  Helpers::SetAllPixelsToValue(projectedImage.GetPointer(), emptyVector);

  itk::ImageRegionConstIterator<PTXImage::XYZImageType> xyzImageIterator(xyzImage, xyzImage->GetLargestPossibleRegion());

  unsigned int badPoints = 0;

  // Iterate over the scan points image and track where each projects in the 'projectedImage'
  while(!xyzImageIterator.IsAtEnd())
    {
    if(ptxImage.GetPTXPixel(xyzImageIterator.GetIndex()).Valid)
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

      if(!colorImage->GetLargestPossibleRegion().IsInside(projectedPixel))
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
      } // end if valid

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
        } // end loop over points projected to this pixel

      // Color all points within a tolerance of the minimum depth
      PTXImage::RGBImageType::PixelType color;
      color = colorImage->GetPixel(projectedImageIterator.GetIndex());
      for(unsigned int i = 0; i < projectedImageIterator.Get().size(); ++i)
        {
        itk::Index<2> currentPixel = projectedImageIterator.Get()[i];
        PTXPixel ptxPixel = ptxImage.GetPTXPixel(currentPixel);
        if(ptxPixel.Valid)
          {
          if(fabs(ptxPixel.GetDepth() - minDepth) < .1)
            {
            colorImage->SetPixel(currentPixel, color);
            validityMask->SetPixel(currentPixel, 255);
            } // end if in tolerance
          } // end if valid
        } // end for

      } // end else over > 0 projections

    ++projectedImageIterator;
    } // end while over whole image


  typedef  itk::ImageFileWriter< PTXImage::MaskImageType > MaskWriterType;
  MaskWriterType::Pointer maskWriter = MaskWriterType::New();
  maskWriter->SetFileName("ValidColorMask.png");
  maskWriter->SetInput(validityMask);
  maskWriter->Update();

  ptxImage.ReplaceValidity(validityMask);

  ptxImage.ReplaceRGB(colorImage);

  std::cout << "Output has " << ptxImage.CountValidPoints() << " valid points." << std::endl;

}

/*
void ResectionNaive()
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
  PTXImage ptxImage = PTXReader::Read(ptxFileName);

  PTXImage::XYZImageType::Pointer xyzImage = ptxImage.GetXYZImage();

  PTXImage::RGBImageType::Pointer colorImage = PTXImage::RGBImageType::New();
  ptxImage.CreateRGBImage(colorImage);
  //std::cout << "colorImage: " << colorImage->GetLargestPossibleRegion() << std::endl;

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
  unsigned int badPoints = 0;
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
      //std::cout << "Point does not project to image!" << std::endl;
      badPoints++;
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

  std::cout << "There were " << badPoints << " points that did not project to inside of the image!" << std::endl;
  ptxImage.ReplaceRGB(colorImage);

  //FilePrefix vtpPrefix("outputPointCloud");
  //ptxImage.WritePointCloud(vtpPrefix);

  FilePrefix ptxPrefix(outputFileName);
  ptxImage.WritePTX(ptxPrefix);

}*/
