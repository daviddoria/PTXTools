#include "PTXImage.h"
#include "Frame.h"

// ITK
#include "itkCovariantVector.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkPoint.h"
#include "itkAzimuthElevationToCartesianTransform.h"
#include "itkShrinkImageFilter.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkConstantBoundaryCondition.h"

// VTK
#include <vtkAppendPolyData.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkPlaneSource.h>
#include <vtkLineSource.h>
#include <vtkPlane.h>
#include <vtkMath.h>
#include <vtkUnsignedCharArray.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkXMLPolyDataWriter.h>

// STL
#include <fstream>
#include <sstream>
#include <string>

PTXImage::PTXImage()
{
  // Create the main image
  this->FullImage = FullImageType::New();
}

void PTXImage::CountInvalidPoints()
{
  itk::ImageRegionConstIterator<FullImageType> fullImageIterator(this->FullImage, this->FullImage->GetLargestPossibleRegion());

  unsigned int numberOfInvalidPoints = 0;
  while(!fullImageIterator.IsAtEnd())
    {
    PTXPixel fullPixel = fullImageIterator.Get();

    if(!fullPixel.Valid)
      {
      numberOfInvalidPoints++;
      }

    ++fullImageIterator;
    }

  std::cout << "There are " << numberOfInvalidPoints << " invalid points." << std::endl;
}

void PTXImage::WriteDepthThresholdMask(std::string& filename, float depthThreshold)
{
  MaskImageType::Pointer mask = MaskImageType::New();
  mask->SetRegions(this->FullImage->GetLargestPossibleRegion());
  mask->Allocate();
  mask->FillBuffer(0);

  // Setup iterators
  itk::ImageRegionIterator<MaskImageType> maskIterator(mask, mask->GetLargestPossibleRegion());

  itk::ImageRegionConstIterator<FullImageType> fullImageIterator(this->FullImage, this->FullImage->GetLargestPossibleRegion());

  while(!fullImageIterator.IsAtEnd())
    {
    float depth = fullImageIterator.Get().GetDepth();

    if(depth > depthThreshold)
      {
      maskIterator.Set(0); // zero (black) means "valid"
      }
    else
      {
      maskIterator.Set(255); // nonzero means "unknown/invalid"
      }

    ++fullImageIterator;
    ++maskIterator;
    }

  std::stringstream ss;
  ss << filename << "_DepthThresholdMask.png";

  typedef  itk::ImageFileWriter<MaskImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(ss.str());
  writer->SetInput(mask);
  writer->Update();
}

void PTXImage::WriteInvalidMask(std::string& filename)
{
  MaskImageType::Pointer mask = MaskImageType::New();
  mask->SetRegions(this->FullImage->GetLargestPossibleRegion());
  mask->Allocate();
  mask->FillBuffer(0);

  // Setup iterators
  itk::ImageRegionIterator<MaskImageType> maskIterator(mask, mask->GetLargestPossibleRegion());

  itk::ImageRegionConstIterator<FullImageType> fullImageIterator(this->FullImage, this->FullImage->GetLargestPossibleRegion());

  while(!fullImageIterator.IsAtEnd())
    {
    PTXPixel fullPixel = fullImageIterator.Get();

    if(fullPixel.Valid)
      {
      maskIterator.Set(0); // zero (black) means "valid"
      }
    else
      {
      maskIterator.Set(255); // nonzero means "unknown/invalid"
      }

    ++fullImageIterator;
    ++maskIterator;
    }

  std::stringstream ss;
  ss << filename << "_ValidityMask.png";

  typedef  itk::ImageFileWriter<MaskImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(ss.str());
  writer->SetInput(mask);
  writer->Update();
}

void PTXImage::ReadFile(std::string filename)
{
  // Attempt to open the file
  std::ifstream infile;
  infile.open(filename.c_str());

  // Verify that the file was opened correctly
  if(!infile)
    {
    std::cout << "Could not open file " << filename << "!" << std::endl;
    return;
    }

  // Read the header
  std::string line;

  unsigned int numberOfThetaPoints; // aka numberOfColumns
  unsigned int numberOfPhiPoints; // aka numberOfRows

  getline(infile, line);
  std::stringstream(line) >> numberOfThetaPoints;

  getline(infile, line);
  std::stringstream(line) >> numberOfPhiPoints;

  std::cout << "PhiPoints: " << numberOfPhiPoints << std::endl;
  std::cout << "ThetaPoints: " << numberOfThetaPoints << std::endl;

  // Setup the image
  itk::Size<2> size;
  size[0] = numberOfThetaPoints;
  size[1] = numberOfPhiPoints;

  itk::Index<2> start;
  start.Fill(0);

  itk::ImageRegion<2> region(start, size);

  this->FullImage->SetRegions(region);
  this->FullImage->Allocate();

  // Skip 8 lines (identity matrices)
  for(int i = 0; i < 8; i++)
    {
    getline(infile, line);
    }

  for(unsigned int theta = 0; theta < numberOfThetaPoints; theta++)
    {
    for(unsigned int phi = 0; phi < numberOfPhiPoints; phi++)
      {
      itk::Index<2> pixelIndex;
      pixelIndex[0] = theta;
      pixelIndex[1] = phi;

      getline(infile, line);

      float coordinate[3];
      float intensity;
      //unsigned char color[3];
      int color[3]; // must read these from the string stream as an int

      std::stringstream parsedLine(line);
      parsedLine >> coordinate[0] >> coordinate[1] >> coordinate[2] >> intensity
                 >> color[0] >> color[1] >> color[2];
      //std::cout << color[0] << " " << color[1] << " " << color[2] << std::endl;

      // Create the pixel with all of the properties
      PTXPixel pixel;
      pixel.X = coordinate[0];
      pixel.Y = coordinate[1];
      pixel.Z = coordinate[2];
      pixel.Intensity = intensity;
      pixel.R = color[0];
      pixel.G = color[1];
      pixel.B = color[2];

      // Check for a particular value (0 0 0 0.5 0 0 0) which indicates that the scanner did not receive a valid return.
      if(pixel.X == 0 && pixel.Y == 0 && pixel.Z == 0 &&
        intensity == 0.50 && pixel.R == 0 && pixel.G == 0 && pixel.B == 0)
        {
        pixel.Valid = false;
        }
      else
        {
        pixel.Valid = true;
        }

      // Set the pixel in the image
      this->FullImage->SetPixel(pixelIndex, pixel);
      //std::cout << "Pixel " << pixelIndex << " Valid? " << pixel.Valid << std::endl;
      }//end phi for loop
    }// end theta for loop

  // Close the input file
  infile.close();

}

void PTXImage::CreateRGBImage(RGBImageType::Pointer image)
{
  // Setup the image
  image->SetRegions(this->FullImage->GetLargestPossibleRegion());
  image->Allocate();
  image->FillBuffer(itk::NumericTraits<RGBImageType::PixelType>::Zero);

  // Setup iterators
  itk::ImageRegionIterator<RGBImageType> rgbImageIterator(image, image->GetLargestPossibleRegion());
  rgbImageIterator.GoToBegin();

  itk::ImageRegionIterator<FullImageType> fullImageIterator(this->FullImage, this->FullImage->GetLargestPossibleRegion());
  fullImageIterator.GoToBegin();

  // Traverse the full image, extracting the color information and saving it in the RGB image
  while(!rgbImageIterator.IsAtEnd())
    {
    PTXPixel fullPixel = fullImageIterator.Get();

    itk::CovariantVector<unsigned char, 3> rgbPixel;

    if(fullPixel.Valid)
      {
      rgbPixel[0] = fullPixel.R;
      rgbPixel[1] = fullPixel.G;
      rgbPixel[2] = fullPixel.B;
      }
    else
      {
      // Set the pixel to bright green
      rgbPixel[0] = 0;
      //rgbPixel[1] = 255;
      rgbPixel[1] = 0;
      rgbPixel[2] = 0;
      }

    rgbImageIterator.Set(rgbPixel);

    ++rgbImageIterator;
    ++fullImageIterator;
    }

}

void PTXImage::WriteRGBImage(std::string filename)
{
  // This a convenience function which simply calls CreateRGBImage and then writes the result to a file
  RGBImageType::Pointer image = RGBImageType::New();
  CreateRGBImage(image);

  std::stringstream ss;
  ss << filename << "_RGB.png";

  typedef  itk::ImageFileWriter< RGBImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(ss.str());
  writer->SetInput(image);
  writer->Update();

}

void PTXImage::CreatePointCloud(vtkSmartPointer<vtkPolyData> pointCloud)
{
  // Create point and color arrays
  vtkSmartPointer<vtkPoints> points =
    vtkSmartPointer<vtkPoints>::New();

  vtkSmartPointer<vtkUnsignedCharArray> colors =
    vtkSmartPointer<vtkUnsignedCharArray>::New();
  colors->SetNumberOfComponents(3);
  colors->SetName("Colors");

  // Iterate through the full image extracting the coordinate and color information and adding them to their respective arrays
  itk::ImageRegionConstIteratorWithIndex<FullImageType> imageIterator(this->FullImage, this->FullImage->GetLargestPossibleRegion());

  while(!imageIterator.IsAtEnd())
    {
    PTXPixel pixel = imageIterator.Get();
    unsigned char rgb[3];
    rgb[0] = pixel.R;
    rgb[1] = pixel.G;
    rgb[2] = pixel.B;

    if(pixel.Valid)
      {
      colors->InsertNextTupleValue(rgb);
      points->InsertNextPoint(pixel.X, pixel.Y, pixel.Z);
      }

    ++imageIterator;
    }

  // Combine the coordinates and colors into a PolyData
  vtkSmartPointer<vtkPolyData> polydata =
    vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints(points);
  polydata->GetPointData()->SetScalars(colors);

  // Create a vertex at each point
  vtkSmartPointer<vtkVertexGlyphFilter> vertexGlyphFilter =
    vtkSmartPointer<vtkVertexGlyphFilter>::New();
  vertexGlyphFilter->AddInput(polydata);
  vertexGlyphFilter->Update();

  pointCloud->ShallowCopy(vertexGlyphFilter->GetOutput());
}

void PTXImage::WritePointCloud(std::string filename)
{
  // This a convenience function which simply calls CreatePointCloud and then writes the result to a file
  std::stringstream ss;
  ss << filename << ".vtp";

  vtkSmartPointer<vtkPolyData> pointCloud = vtkSmartPointer<vtkPolyData>::New();
  CreatePointCloud(pointCloud);

  vtkSmartPointer<vtkXMLPolyDataWriter> writer =
    vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName(ss.str().c_str());
  writer->SetInputConnection(pointCloud->GetProducerPort());
  writer->Write();
}

void PTXImage::CreateDepthImage(FloatImageType::Pointer image)
{
  // Setup the image
  image->SetRegions(this->FullImage->GetLargestPossibleRegion());
  image->Allocate();
  image->FillBuffer(0);

  // Setup iterators
  itk::ImageRegionConstIteratorWithIndex<FullImageType> imageIterator(this->FullImage, this->FullImage->GetLargestPossibleRegion());
  itk::ImageRegionIterator<FloatImageType> depthImageIterator(image, image->GetLargestPossibleRegion());

  // Traverse the full image extracting the depth information and saving it in the depth image
  while(!imageIterator.IsAtEnd())
    {
    depthImageIterator.Set(imageIterator.Get().GetDepth());

    ++imageIterator;
    ++depthImageIterator;
    }

}

void PTXImage::WriteDepthImage(std::string filePrefix)
{
  // This a convenience function which simply calls CreateDepthImage and then writes the result to png (scaled) and mhd (unscaled) files.
  // 'filePrefix' is, for example "filename" which will be used to internally produce the filenames "filename.mhd" and "filename.png"

  FloatImageType::Pointer image = FloatImageType::New();
  CreateDepthImage(image);

  std::stringstream ssMHD;
  ssMHD << filePrefix << "Depth.mhd";
  typedef  itk::ImageFileWriter<FloatImageType> MHDWriterType;
  MHDWriterType::Pointer mhdWriter = MHDWriterType::New();
  mhdWriter->SetFileName(ssMHD.str());
  mhdWriter->SetInput(image);
  mhdWriter->Update();

  typedef itk::Image<unsigned char, 2> UnsignedCharImageType;
  typedef itk::RescaleIntensityImageFilter<FloatImageType, UnsignedCharImageType > RescaleFilterType;
  RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
  rescaleFilter->SetInput(image);
  rescaleFilter->SetOutputMinimum(0);
  rescaleFilter->SetOutputMaximum(255);
  rescaleFilter->Update();

  std::stringstream ssPNG;
  ssPNG << filePrefix << "Depth.png";

  typedef  itk::ImageFileWriter< UnsignedCharImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(ssPNG.str());
  writer->SetInput(rescaleFilter->GetOutput());
  writer->Update();

}


void PTXImage::CreateIntensityImage(FloatImageType::Pointer image)
{
  // Setup the image
  image->SetRegions(this->FullImage->GetLargestPossibleRegion());
  image->Allocate();
  image->FillBuffer(0);

  // Setup iterators
  itk::ImageRegionConstIteratorWithIndex<FullImageType> imageIterator(this->FullImage, this->FullImage->GetLargestPossibleRegion());
  itk::ImageRegionIterator<FloatImageType> intensityImageIterator(image, image->GetLargestPossibleRegion());

  // Traverse the full image extracting the depth information and saving it in the depth image
  while(!imageIterator.IsAtEnd())
    {
    PTXPixel pixel = imageIterator.Get();
    if(pixel.Valid)
      {
      intensityImageIterator.Set(pixel.Intensity);
      }
    else
      {
      intensityImageIterator.Set(0);
      }

    ++imageIterator;
    ++intensityImageIterator;
    }

}

void PTXImage::WriteIntensityImage(std::string filePrefix)
{
  // This a convenience function which simply calls CreateIntensityImage and then writes the result to png (scaled) and mhd (unscaled) files.
  // 'filePrefix' is, for example "filename" which will be used to internally produce the filenames "filename.mhd" and "filename.png"

  FloatImageType::Pointer image = FloatImageType::New();
  CreateIntensityImage(image);

  std::stringstream ssMHD;
  ssMHD << filePrefix << "Intensity.mhd";
  typedef  itk::ImageFileWriter< FloatImageType > MHDWriterType;
  MHDWriterType::Pointer mhdWriter = MHDWriterType::New();
  mhdWriter->SetFileName(ssMHD.str());
  mhdWriter->SetInput(image);
  mhdWriter->Update();

  typedef itk::Image<unsigned char, 2> UnsignedCharImageType;
  typedef itk::RescaleIntensityImageFilter< FloatImageType, UnsignedCharImageType > RescaleFilterType;
  RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
  rescaleFilter->SetInput(image);
  rescaleFilter->SetOutputMinimum(0);
  rescaleFilter->SetOutputMaximum(255);
  rescaleFilter->Update();

  std::stringstream ssPNG;
  ssPNG << filePrefix << "Intensity.png";

  typedef  itk::ImageFileWriter< UnsignedCharImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(ssPNG.str());
  writer->SetInput(rescaleFilter->GetOutput());
  writer->Update();

}

void PTXImage::ReplaceDepth(itk::Image<float, 2>::Pointer depthImage)
{
  // This function allows the depth map to be modified externally and the new map applied to the grid

  CountInvalidPoints();

  // Setup iterators
  itk::ImageRegionIterator<FullImageType> imageIterator(this->FullImage, this->FullImage->GetLargestPossibleRegion());
  itk::ImageRegionConstIterator<itk::Image<float, 2> > newDepthIterator(depthImage, depthImage->GetLargestPossibleRegion());

  itk::Point<float, 3> origin;
  origin.Fill(0);

  while(!imageIterator.IsAtEnd())
    {
    // Get the old point
    PTXPixel pixel = imageIterator.Get();
    itk::Point<float, 3> oldPoint;
    oldPoint[0] = pixel.X;
    oldPoint[1] = pixel.Y;
    oldPoint[2] = pixel.Z;

    // For testing only, compute the old depth
    double oldDepth = origin.EuclideanDistanceTo(oldPoint);
    std::cout << "Old depth: " << oldDepth << " New depth: " << newDepthIterator.Get() << std::endl;

    // Get the vector from the origin (scanner location) and the old point
    itk::Vector<float, 3> unitVector = oldPoint - origin;

    // Get a unit vector in the direction of the old point
    unitVector.Normalize();

    // Compute the new point from the vector and the new depth
    itk::Point<float, 3> newPoint = origin + unitVector * newDepthIterator.Get();

    std::cout << "Old point: " << oldPoint << " New point: " << newPoint << std::endl;
    // Save the new point in the PTXPixel
    pixel.X = newPoint[0];
    pixel.Y = newPoint[1];
    pixel.Z = newPoint[2];

    imageIterator.Set(pixel);

    ++imageIterator;
    ++newDepthIterator;
    }

}


void PTXImage::ReplaceRGBD(itk::Image<itk::CovariantVector<float, 4>, 2>::Pointer rgbd)
{
  // This function allows the color and depth to be modified externally and the new map applied to the grid

  CountInvalidPoints();

  ComputeAverageDeltaPhi();
  ComputeAverageDeltaTheta();

  // Setup iterators
  itk::ImageRegionIterator<FullImageType> imageIterator(this->FullImage, this->FullImage->GetLargestPossibleRegion());
  itk::ImageRegionConstIterator<itk::Image<itk::CovariantVector<float, 4>, 2> > rgbdIterator(rgbd, rgbd->GetLargestPossibleRegion());

  itk::Point<float, 3> origin;
  origin.Fill(0);

  unsigned int madeValid = 0;

  while(!imageIterator.IsAtEnd())
    {
    // Get the old point
    PTXPixel pixel = imageIterator.Get();

    // Copy the color from the RGBD image
    pixel.R = rgbdIterator.Get()[0];
    pixel.G = rgbdIterator.Get()[1];
    pixel.B = rgbdIterator.Get()[2];

    itk::Point<float, 3> oldPoint;

    // If the point was invalid, construct a vector in its direction incase the depth is now valid
    bool previousValidity = pixel.Valid;

    if(pixel.Valid == false)
      {
      madeValid++;
      pixel.Valid = true;

      oldPoint = ApproximateOldPoint(imageIterator.GetIndex());
      std::cout << "new oldPoint: " << oldPoint << std::endl;

      // For debugging, make new pixels bright green
//       pixel.R = 0;
//       pixel.G = 255;
//       pixel.B = 0;

      // For debugging, skip invalid pixels
      ++imageIterator;
      ++rgbdIterator;
      continue;
      }
    else
      {
      // For debugging, turn off the pixels which are not newly valid
      //pixel.Valid = false; //!!!

      oldPoint[0] = pixel.X;
      oldPoint[1] = pixel.Y;
      oldPoint[2] = pixel.Z;
      }

    // Get the vector from the origin (scanner location) and the old point
    itk::Vector<float, 3> unitVector = oldPoint - origin;

    // Get a unit vector in the direction of the old point
    unitVector.Normalize();

    // Compute the new point from the vector and the new depth
    float multiplier = 1.0;
    if(!previousValidity)
      {
      multiplier = -1.0;
      }
    float depth = multiplier * rgbdIterator.Get()[3]; // This works, but why is the negative necessary?
    itk::Point<float, 3> newPoint = origin + unitVector * depth;

    // Save the new point in the PTXPixel
    pixel.X = newPoint[0];
    pixel.Y = newPoint[1];
    pixel.Z = newPoint[2];

    imageIterator.Set(pixel);
    //std::cout << pixel << std::endl;

    ++imageIterator;
    ++rgbdIterator;
    }

  std::cout << madeValid << " pixels were made valid." << std::endl;
  CountInvalidPoints();
}

itk::Point<float, 3> PTXImage::ApproximateOldPoint(itk::Index<2> pixel)
{
  // This function creates a point unit distance from the origin in the approximate direction the old point would have been aquired.
  float phi = ApproximatePhi(pixel);
  float theta = ApproximateTheta(pixel);

  typedef itk::Point<float, 3> PointType;
  PointType azEl;
  azEl[0] = theta;
  azEl[1] = phi;
  azEl[2] = 1;

  typedef itk::AzimuthElevationToCartesianTransform< float, 3 >
    AzimuthElevationToCartesian;
  AzimuthElevationToCartesian::Pointer azimuthElevation =
    AzimuthElevationToCartesian::New();
  return azimuthElevation->TransformAzElToCartesian(azEl);
}

void PTXImage::ApplyMask(itk::Image<unsigned char, 2>::Pointer mask)
{
  // Blank the PTX image in areas where mask is non-zero
  itk::ImageRegionIterator<FullImageType> imageIterator(this->FullImage, this->FullImage->GetLargestPossibleRegion());
  itk::ImageRegionConstIterator<itk::Image<unsigned char, 2> > maskIterator(mask, mask->GetLargestPossibleRegion());

  while(!imageIterator.IsAtEnd())
    {
    if(maskIterator.Get())
      {
      PTXPixel pixel = imageIterator.Get();
      pixel.Valid = false;
      imageIterator.Set(pixel);
      }
    ++imageIterator;
    ++maskIterator;
    }
}

void PTXImage::CreateRGBDIImage(RGBDIImageType::Pointer image)
{
  // Create the 5 channels
  FloatImageType::Pointer intensityImage = FloatImageType::New();
  CreateIntensityImage(intensityImage);

  FloatImageType::Pointer depthImage = FloatImageType::New();
  CreateDepthImage(depthImage);

  RGBImageType::Pointer rgbImage = RGBImageType::New();
  CreateRGBImage(rgbImage);

  // Setup the image
  image->SetRegions(this->FullImage->GetLargestPossibleRegion());
  image->Allocate();
  image->FillBuffer(itk::NumericTraits<RGBDIImageType::PixelType>::Zero);

  // Setup iterators
  itk::ImageRegionIterator<RGBDIImageType> imageIterator(image, image->GetLargestPossibleRegion());
  itk::ImageRegionConstIterator<FloatImageType> intensityImageIterator(intensityImage, intensityImage->GetLargestPossibleRegion());
  itk::ImageRegionConstIterator<FloatImageType> depthImageIterator(depthImage, depthImage->GetLargestPossibleRegion());
  itk::ImageRegionConstIterator<RGBImageType> rgbImageIterator(rgbImage, rgbImage->GetLargestPossibleRegion());

  // Copy the images into their respective channels
  while(!imageIterator.IsAtEnd())
    {

    RGBDIImageType::PixelType pixel;
    pixel[0] = rgbImageIterator.Get()[0];
    pixel[1] = rgbImageIterator.Get()[1];
    pixel[2] = rgbImageIterator.Get()[2];
    pixel[3] = depthImageIterator.Get();
    pixel[4] = intensityImageIterator.Get();

    imageIterator.Set(pixel);

    ++imageIterator;
    ++intensityImageIterator;
    ++depthImageIterator;
    ++rgbImageIterator;
    }
}

void PTXImage::WriteRGBDIImage(std::string filePrefix)
{
  RGBDIImageType::Pointer image = RGBDIImageType::New();
  CreateRGBDIImage(image);

  std::stringstream ss;
  ss << filePrefix << "_RGBDI.mhd";

  typedef  itk::ImageFileWriter< RGBDIImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(ss.str());
  writer->SetInput(image);
  writer->Update();
}

void PTXImage::CreateRGBDImage(RGBDImageType::Pointer image)
{
  // Create the 4 channels
  FloatImageType::Pointer depthImage = FloatImageType::New();
  CreateDepthImage(depthImage);

  RGBImageType::Pointer rgbImage = RGBImageType::New();
  CreateRGBImage(rgbImage);

  // Setup the image
  image->SetRegions(this->FullImage->GetLargestPossibleRegion());
  image->Allocate();
  image->FillBuffer(itk::NumericTraits<RGBDImageType::PixelType>::Zero);

  // Setup iterators
  itk::ImageRegionIterator<RGBDImageType> imageIterator(image, image->GetLargestPossibleRegion());
  itk::ImageRegionConstIterator<FloatImageType> depthImageIterator(depthImage, depthImage->GetLargestPossibleRegion());
  itk::ImageRegionConstIterator<RGBImageType> rgbImageIterator(rgbImage, rgbImage->GetLargestPossibleRegion());

  // Copy the images into their respective channels
  while(!imageIterator.IsAtEnd())
    {

    RGBDImageType::PixelType pixel;
    pixel[0] = rgbImageIterator.Get()[0];
    pixel[1] = rgbImageIterator.Get()[1];
    pixel[2] = rgbImageIterator.Get()[2];
    pixel[3] = depthImageIterator.Get();

    imageIterator.Set(pixel);

    ++imageIterator;
    ++depthImageIterator;
    ++rgbImageIterator;
    }
}

void PTXImage::WriteRGBDImage(std::string filePrefix)
{
  RGBDImageType::Pointer image = RGBDImageType::New();
  CreateRGBDImage(image);

  std::stringstream ss;
  ss << filePrefix << "_RGBD.mhd";

  typedef  itk::ImageFileWriter< RGBDImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(ss.str());
  writer->SetInput(image);
  writer->Update();
}

void PTXImage::ComputeAverageDeltaTheta()
{
  itk::ImageRegion<2> region = this->FullImage->GetLargestPossibleRegion();

  std::vector<float> deltaThetas;

  itk::Index<2> lastIndex;
  lastIndex[0] = -1;
  lastIndex[1] = -1;

  for(unsigned int row = 0; row < region.GetSize()[1]; row++)
    {
    for(unsigned int col = 0; col < region.GetSize()[0]; col++)
      {
      itk::Index<2> index;
      index[0] = col;
      index[1] = row;

      if(!this->FullImage->GetPixel(index).Valid || // the current pixel is invalid
        !this->FullImage->GetPixel(lastIndex).Valid || // the previous pixel was invalid
        lastIndex[1] != index[1] // if we are not on the same row as the last pixel
        )
        {
        // do nothing
        }
      else
        {
        deltaThetas.push_back(GetTheta(index) - GetTheta(lastIndex));
        }

      lastIndex = index;
      }
    }

  float total = 0.0;
  for(unsigned int i = 0; i < deltaThetas.size(); i++)
    {
    total += deltaThetas[i];
    }

  if(deltaThetas.size() == 0)
    {
    std::cerr << "There are no valid theta neighbors!" << std::endl;
    exit(-1);
    }

  this->AverageDeltaTheta = total/static_cast<float>(deltaThetas.size());

  std::cout << "AverageDeltaTheta: " << this->AverageDeltaTheta << std::endl;
}

void PTXImage::ComputeAverageDeltaPhi()
{
  itk::ImageRegion<2> region = this->FullImage->GetLargestPossibleRegion();

  std::vector<float> deltaPhis;

  itk::Index<2> lastIndex;
  lastIndex[0] = -1;
  lastIndex[1] = -1;

  for(unsigned int col = 0; col < region.GetSize()[0]; col++)
    {
    for(unsigned int row = 0; row < region.GetSize()[1]; row++)
      {
      itk::Index<2> index;
      index[0] = col;
      index[1] = row;

      if(!this->FullImage->GetPixel(index).Valid || // the current pixel is invalid
        !this->FullImage->GetPixel(lastIndex).Valid || // the previous pixel is invalid
        lastIndex[0] != index[0] // if we are not on the same column as the last pixel
        )
        {
        // do nothing
        }
      else
        {
        deltaPhis.push_back(GetPhi(index) - GetPhi(lastIndex));
        }

      lastIndex = index;
      }
    }

  float total = 0.0;
  for(unsigned int i = 0; i < deltaPhis.size(); i++)
    {
    total += deltaPhis[i];
    }

  if(deltaPhis.size() == 0)
    {
    std::cerr << "There are no valid phi neighbors!" << std::endl;
    exit(-1);
    }

  this->AverageDeltaPhi = total/static_cast<float>(deltaPhis.size());

  std::cout << "AverageDeltaPhi: " << this->AverageDeltaPhi << std::endl;
}

float PTXImage::GetPhi(itk::Index<2> index)
{
  if(!this->FullImage->GetPixel(index).Valid)
    {
    std::cerr << "Cannot GetPhi on an invalid pixel!" << std::endl;
    exit(-1);
    }

  typedef itk::Point<float, 3> PointType;
  PointType cartesian;
  cartesian[0] = this->FullImage->GetPixel(index).X;
  cartesian[1] = this->FullImage->GetPixel(index).Y;
  cartesian[2] = this->FullImage->GetPixel(index).Z;

  typedef itk::AzimuthElevationToCartesianTransform< float, 3 >
    AzimuthElevationToCartesian;
  AzimuthElevationToCartesian::Pointer azimuthElevation =
    AzimuthElevationToCartesian::New();

  return azimuthElevation->TransformCartesianToAzEl(cartesian)[1];
}

float PTXImage::GetTheta(itk::Index<2> index)
{
  if(!this->FullImage->GetPixel(index).Valid)
    {
    std::cerr << "Cannot GetPhi on an invalid pixel!" << std::endl;
    exit(-1);
    }

  typedef itk::Point<float, 3> PointType;
  PointType cartesian;
  cartesian[0] = this->FullImage->GetPixel(index).X;
  cartesian[1] = this->FullImage->GetPixel(index).Y;
  cartesian[2] = this->FullImage->GetPixel(index).Z;

  typedef itk::AzimuthElevationToCartesianTransform< float, 3 >
    AzimuthElevationToCartesian;
  AzimuthElevationToCartesian::Pointer azimuthElevation =
    AzimuthElevationToCartesian::New();

  return azimuthElevation->TransformCartesianToAzEl(cartesian)[0];
}

float PTXImage::ApproximateTheta(itk::Index<2> pixel)
{
  itk::Offset<2> offset;
  offset[0] = 1;
  offset[1] = 0;

  itk::Index<2> nearestPixel = FindNearestValidPixel(pixel, offset);

  float theta = GetTheta(nearestPixel);

  return theta + (pixel[0] - nearestPixel[0])*this->AverageDeltaTheta;
}

float PTXImage::ApproximatePhi(itk::Index<2> pixel)
{
  itk::Offset<2> offset;
  offset[0] = 0;
  offset[1] = 1;

  itk::Index<2> nearestPixel = FindNearestValidPixel(pixel, offset);

  float phi = GetPhi(nearestPixel);

  return phi + (pixel[1] - nearestPixel[1])*this->AverageDeltaPhi;
}

itk::Index<2> PTXImage::FindNearestValidPixel(itk::Index<2> pixel, itk::Offset<2> offset)
{
  // This function finds the nearest valid pixel along a row or column (specified by offset) of an image
  itk::Index<2> currentPixel = pixel;

  itk::Size<2> size = this->FullImage->GetLargestPossibleRegion().GetSize();

  // Step forward
  unsigned int pixelCounter = 0;
  while(this->FullImage->GetLargestPossibleRegion().IsInside(currentPixel+offset))
  {
    currentPixel += offset;
    //std::cout << "currentPixel: " << currentPixel << " Valid? " << this->FullImage->GetPixel(currentPixel).Valid << std::endl;
    pixelCounter++;
    if(this->FullImage->GetPixel(currentPixel).Valid)
      {
      std::cout << "Closest valid pixel to " << pixel << " is " << currentPixel << std::endl;
      return currentPixel;
      }
  }

  std::cout << "Searched " << pixelCounter << " pixels." << std::endl;

  std::cout << "Offset was " << offset << std::endl;
  offset[0] *= -1;
  offset[1] *= -1;
  std::cout << "Offset is now " << offset << std::endl;

  currentPixel = pixel; // reset to function input

  // Step backward
  // This is EXACTLY the same loop, but the offset has switched directions
  while(this->FullImage->GetLargestPossibleRegion().IsInside(currentPixel+offset))
  {
    currentPixel += offset;
    //std::cout << "currentPixel: " << currentPixel << " Valid? " << this->FullImage->GetPixel(currentPixel).Valid << std::endl;
    pixelCounter++;
    if(this->FullImage->GetPixel(currentPixel).Valid)
      {
      std::cout << "Closest valid pixel to " << pixel << " is " << currentPixel << std::endl;
      return currentPixel;
      }
  }

  //std::cout << "Searched " << pixelCounter << " pixels." << std::endl;

  // If we get here, an entire row or column did not have a valid pixel. We cannot deal with images this poor.
  std::cerr << "No neighbor was found for " << pixel << " in the " << offset << " direction. This should never happen!" << std::endl;
  exit(-1);
  return currentPixel; // So the compiler doesn't complain that all paths don't return a value
}

PTXImage PTXImage::Downsample(const unsigned int factor)
{
  typedef itk::ShrinkImageFilter <FullImageType, FullImageType> ShrinkImageFilterType;
  ShrinkImageFilterType::Pointer shrinkFilter = ShrinkImageFilterType::New();
  shrinkFilter->SetInput(this->FullImage);
  shrinkFilter->SetShrinkFactor(0, factor);
  shrinkFilter->SetShrinkFactor(1, factor);
  shrinkFilter->Update();

  PTXImage output;
  DeepCopy<FullImageType>(shrinkFilter->GetOutput(), output.FullImage);

  return output;
}

void PTXImage::WritePTX(const std::string filename)
{
  std::ofstream fout(filename.c_str());

  unsigned int numberOfThetaPoints = this->FullImage->GetLargestPossibleRegion().GetSize()[0];
  unsigned int numberOfPhiPoints = this->FullImage->GetLargestPossibleRegion().GetSize()[1];

  fout <<  numberOfThetaPoints << std::endl;
  fout << numberOfPhiPoints << std::endl;
  fout << "0 0 0" << std::endl
       << "1 0 0" << std::endl
       << "0 1 0" << std::endl
       << "0 0 1" << std::endl
       << "1 0 0 0" << std::endl
       << "0 1 0 0" << std::endl
       << "0 0 1 0" << std::endl
       << "0 0 0 1" << std::endl;

  for(unsigned int theta = 0; theta < numberOfThetaPoints; theta++)
    {
    for(unsigned int phi = 0; phi < numberOfPhiPoints; phi++)
      {
      itk::Index<2> pixelIndex;
      pixelIndex[0] = theta;
      pixelIndex[1] = phi;

      PTXPixel pixel = this->FullImage->GetPixel(pixelIndex);
      fout << pixel.X << " " << pixel.Y << " " << pixel.Z << " " << pixel.Intensity
           << " " << static_cast<int>(pixel.R) << " " << static_cast<int>(pixel.G)
           << " " << static_cast<int>(pixel.B) << std::endl;
      }
    }

  fout.close();
}

void PTXImage::ComputeWeightedDepthLaplacian(const std::string filename)
{
  typedef itk::Image<float, 2> FloatImageType;
  FloatImageType::Pointer output = FloatImageType::New();
  output->SetRegions(this->FullImage->GetLargestPossibleRegion());
  output->Allocate();

  itk::Size<2> radius;
  radius.Fill(1);
  itk::ConstNeighborhoodIterator<FullImageType> imageIterator(radius, this->FullImage, this->FullImage->GetLargestPossibleRegion());

  FloatImageType::Pointer depthImage = FloatImageType::New();
  CreateDepthImage(depthImage);

  itk::ImageRegionIterator<FloatImageType> outputIterator(output, output->GetLargestPossibleRegion());

  while(!imageIterator.IsAtEnd())
    {
    if(!imageIterator.GetCenterPixel().Valid)
      {
      ++imageIterator;
      ++outputIterator;
      continue;
      }

    std::vector<itk::Offset<2> > offsets;

    itk::Offset<2> left;
    left[0] = -1;
    left[1] = 0;
    offsets.push_back(left);

    itk::Offset<2> right;
    right[0] = 1;
    right[1] = 0;
    offsets.push_back(right);

    itk::Offset<2> top;
    top[0] = 0;
    top[1] = -1;
    offsets.push_back(top);

    itk::Offset<2> bottom;
    bottom[0] = 0;
    bottom[1] = 1;
    offsets.push_back(bottom);

    std::vector<float> inverseDistances;
    std::vector<float> depths;
    for(unsigned int i = 0; i < offsets.size(); i++)
      {
      if(imageIterator.GetPixel(left).Valid)
        {
        float dist = DistanceBetweenPoints(imageIterator.GetCenterPixel(), imageIterator.GetPixel(offsets[i]));
        inverseDistances.push_back(1.0/dist);
        depths.push_back(depthImage->GetPixel(imageIterator.GetIndex() + offsets[i]));
        }
      }

    // Normalize inverseDistances
    float sum = 0.0;
    for(unsigned int i = 0; i < inverseDistances.size(); i++)
      {
      sum += inverseDistances[i];
      }
    for(unsigned int i = 0; i < inverseDistances.size(); i++)
      {
      inverseDistances[i] /= sum;
      }
    for(unsigned int i = 0; i < inverseDistances.size(); i++)
      {
      inverseDistances[i] *= static_cast<float>(inverseDistances[i]); // The number of pixels used in the Laplacian computation
      }

    // Compute weighted Laplacian
    float total = 0.0;
    total += depths.size(); // The value of the center pixel in the kernel is the number of neighbors used in the comparison
    for(unsigned int i = 0; i < depths.size(); i++)
      {
      total += inverseDistances[i] * depths[i];
      }

    ++imageIterator;
    ++outputIterator;
    }

  typedef  itk::ImageFileWriter<FloatImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(filename);
  writer->SetInput(output);
  writer->Update();
}

PTXImage::VectorType PTXImage::GetPrincipalAxis()
{
  // Get center pixel
  itk::Index<2> centerIndex;
  centerIndex[0] = this->FullImage->GetLargestPossibleRegion().GetSize()[0] / 2;
  centerIndex[1] = this->FullImage->GetLargestPossibleRegion().GetSize()[1] / 2;

  FullImageType::PixelType centerPixel = this->FullImage->GetPixel(centerIndex);

  // Assuming the scanner is at the origin, the direction is simply the coordinate
  typedef itk::CovariantVector<double, 3> VectorType;
  VectorType v;
  v[0] = centerPixel.X;
  v[1] = centerPixel.Y;
  v[2] = centerPixel.Z;

  v.Normalize();

  return v;
}

float PTXImage::DistanceBetweenPoints(PTXPixel a, PTXPixel b)
{
  if(!a.Valid || !b.Valid)
    {
    std::cerr << "Warning: Comparing distance between pixels which are not both valid." << std::endl;
    return 0.0;
    }

  itk::Point<float, 3> p0;
  p0[0] = a.X;
  p0[1] = a.Y;
  p0[2] = a.Z;

  itk::Point<float, 3> p1;
  p1[0] = b.X;
  p1[1] = b.Y;
  p1[2] = b.Z;

  return p0.EuclideanDistanceTo(p1);
}

void PTXImage::WriteProjectionPlane(const std::string filename)
{
  // As always, we assume the scanner is at the origin

  VectorType principalAxis = GetPrincipalAxis();

  double origin[3] = {0,0,0};
  vtkSmartPointer<vtkPlaneSource> planeSource =
    vtkSmartPointer<vtkPlaneSource>::New();
  planeSource->SetCenter(origin);
  planeSource->SetNormal(principalAxis[0], principalAxis[1], principalAxis[2]);
  //planeSource->SetNormal(principalAxis.GetDataPointer()); // only works if the vector is double
  planeSource->Update();

  vtkSmartPointer<vtkLineSource> lineSource =
    vtkSmartPointer<vtkLineSource>::New();
  lineSource->SetPoint1(origin);
  lineSource->SetPoint2(principalAxis[0], principalAxis[1], principalAxis[2]);
  lineSource->Update();

  vtkSmartPointer<vtkAppendPolyData> appendFilter =
    vtkSmartPointer<vtkAppendPolyData>::New();
  appendFilter->AddInputConnection(planeSource->GetOutputPort());
  appendFilter->AddInputConnection(lineSource->GetOutputPort());
  appendFilter->Update();

  vtkSmartPointer<vtkXMLPolyDataWriter> writer =
    vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName(filename.c_str());
  writer->SetInputConnection(appendFilter->GetOutputPort());
  writer->Write();
}

void PTXImage::OrthogonalProjection(std::string filename)
{
  WriteProjectionPlane("projectionPlane.vtp"); // for debugging

  VectorType principalAxis = GetPrincipalAxis();

  // Create a plane through the origin "aimed" at the scanned points
  vtkSmartPointer<vtkPlane> plane =
    vtkSmartPointer<vtkPlane>::New();
  plane->SetOrigin(0.0, 0.0, 0.0);
  plane->SetNormal(principalAxis[0], principalAxis[1], principalAxis[2]);

  itk::ImageRegionConstIterator<FullImageType> imageIterator(this->FullImage, this->FullImage->GetLargestPossibleRegion());

  vtkSmartPointer<vtkPoints> projectedPoints =
    vtkSmartPointer<vtkPoints>::New();

  /*
  while(!imageIterator.IsAtEnd())
    {
    FullImageType::PixelType pixel = imageIterator.Get();

    double point[3] = {pixel.X, pixel.Y, pixel.Z};
    double projected[3];
    plane->ProjectPoint(point, projected);

    projectedPoints->InsertNextPoint(projected);

    ++imageIterator;
    }
  */

  vtkSmartPointer<vtkPolyData> pointCloud =
    vtkSmartPointer<vtkPolyData>::New();
  CreatePointCloud(pointCloud);

  for(unsigned int i = 0; i < static_cast<unsigned int>(pointCloud->GetNumberOfPoints()); i++)
    {
    double point[3];
    pointCloud->GetPoint(i, point);
    double projected[3];
    plane->ProjectPoint(point, projected);

    projectedPoints->InsertNextPoint(projected);
    }

  // Create a polydata of the projected points
  vtkSmartPointer<vtkPolyData> projectedPointsPolydata =
    vtkSmartPointer<vtkPolyData>::New();
  projectedPointsPolydata->SetPoints(projectedPoints);
  projectedPointsPolydata->GetPointData()->SetScalars(pointCloud->GetPointData()->GetScalars());

  {
  // Output projected points for debugging
  vtkSmartPointer<vtkVertexGlyphFilter> vertexGlyphFilter =
    vtkSmartPointer<vtkVertexGlyphFilter>::New();
  vertexGlyphFilter->AddInput(projectedPointsPolydata);
  vertexGlyphFilter->Update();

  vtkSmartPointer<vtkXMLPolyDataWriter> writer =
    vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName("projectedPoints.vtp");
  writer->SetInputConnection(vertexGlyphFilter->GetOutputPort());
  writer->Write();
  }

  // Compute the projection of the vector pointing to the top-center pixel on the projection plane

  // Get top center pixel
  itk::Index<2> topCenterIndex;
  topCenterIndex[0] = this->FullImage->GetLargestPossibleRegion().GetSize()[0] / 2;
  topCenterIndex[1] = this->FullImage->GetLargestPossibleRegion().GetSize()[1] - 1;

  FullImageType::PixelType topCenterPixel = this->FullImage->GetPixel(topCenterIndex);

  if(!topCenterPixel.Valid)
    {
    std::cerr << "topCenterPixel is not valid!" << std::endl;
    exit(-1);
    }

  {
  // for debugging only
  topCenterPixel.R = 255;
  topCenterPixel.G = 0;
  topCenterPixel.B = 0;
  this->FullImage->SetPixel(topCenterIndex, topCenterPixel);
  WriteRGBImage("topCenterColoredRed.png");
  }


  // Assuming the scanner is at the origin, the direction is simply the coordinate
  typedef itk::CovariantVector<double, 3> VectorType;
  VectorType topCenterVector;
  topCenterVector[0] = topCenterPixel.X;
  topCenterVector[1] = topCenterPixel.Y;
  topCenterVector[2] = topCenterPixel.Z;
  topCenterVector.Normalize();

  double v[3];
  v[0] = topCenterPixel.X;
  v[1] = topCenterPixel.Y;
  v[2] = topCenterPixel.Z;

  double up[3];
  plane->ProjectVector(v, up);

  // We will eventually create an image with orientation described by 'up' and the cross product of 'up' and 'principalAxis'

  double horizontalAxis[3];
  double principalAxisArray[3] = {principalAxis[0], principalAxis[1], principalAxis[2]};
  vtkMath::Cross(principalAxisArray, up, horizontalAxis);


  // Create the standard frame
  float standardFrameOrigin[3] = {0,0,0};
  float standardFrameXDirection[3] = {1,0,0};
  float standardFrameYDirection[3] = {0,1,0};
  float standardFrameZDirection[3] = {0,0,1};
  Frame standardFrame(standardFrameOrigin, standardFrameXDirection, standardFrameYDirection, standardFrameZDirection);
  standardFrame.Write("standardFrame.vtp");

  // Get center pixel
  itk::Index<2> centerIndex;
  centerIndex[0] = this->FullImage->GetLargestPossibleRegion().GetSize()[0] / 2;
  centerIndex[1] = this->FullImage->GetLargestPossibleRegion().GetSize()[1] / 2;

  FullImageType::PixelType centerPixel = this->FullImage->GetPixel(centerIndex);

  {
  // for debugging only
  centerPixel.R = 0;
  centerPixel.G = 255;
  centerPixel.B = 0;
  this->FullImage->SetPixel(centerIndex, centerPixel);
  WriteRGBImage("centerColoredGreen.png");
  }

  if(!centerPixel.Valid)
    {
    std::cerr << "CenterPixel is not valid!" << std::endl;
    exit(-1);
    }
  double centerPixelArray[3] = {centerPixel.X, centerPixel.Y, centerPixel.Z};

  double projectedCenter[3];
  plane->ProjectPoint(centerPixelArray, projectedCenter);

  //float originalFrameOrigin[3] = {projectedCenter[0], projectedCenter[1], projectedCenter[2]};
  float originalFrameOrigin[3] = {0,0,0};
  //float originalFrameXDirection[3] = {horizontalAxis[0], horizontalAxis[1], horizontalAxis[2]};
  float originalFrameXDirection[3] = {-horizontalAxis[0], -horizontalAxis[1], -horizontalAxis[2]};
  float originalFrameYDirection[3] = {up[0], up[1], up[2]};
  float originalFrameZDirection[3] = {principalAxisArray[0], principalAxisArray[1], principalAxisArray[2]};
  Frame originalFrame(originalFrameOrigin, originalFrameXDirection, originalFrameYDirection, originalFrameZDirection);
  originalFrame.Write("originalFrame.vtp");

  vtkSmartPointer<vtkTransform> transform =
    vtkSmartPointer<vtkTransform>::New();
  originalFrame.AlignTo(standardFrame, transform);

  vtkSmartPointer<vtkTransformFilter> transformFilter =
    vtkSmartPointer<vtkTransformFilter>::New();
  transformFilter->SetInputConnection(projectedPointsPolydata->GetProducerPort());
  transformFilter->SetTransform(transform);
  transformFilter->Update();

  {
  // Output projected points for debugging
  vtkSmartPointer<vtkVertexGlyphFilter> vertexGlyphFilter =
    vtkSmartPointer<vtkVertexGlyphFilter>::New();
  vertexGlyphFilter->AddInput(transformFilter->GetOutput());
  vertexGlyphFilter->Update();

  vtkSmartPointer<vtkXMLPolyDataWriter> writer =
    vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName("transformedProjectedPoints.vtp");
  writer->SetInputConnection(vertexGlyphFilter->GetOutputPort());
  writer->Write();
  }

  // Compute the corner of the image
  double bounds[6];
  transformFilter->GetOutput()->GetBounds(bounds);

  vtkSmartPointer<vtkTransform> translation =
    vtkSmartPointer<vtkTransform>::New();
  float xmin = bounds[0];
  float ymin = bounds[2];
  translation->Translate(-xmin, -ymin, 0);

  vtkSmartPointer<vtkTransformFilter> shiftTransformFilter =
    vtkSmartPointer<vtkTransformFilter>::New();
  shiftTransformFilter->SetInputConnection(transformFilter->GetOutputPort());
  shiftTransformFilter->SetTransform(translation);
  shiftTransformFilter->Update();

  {
  // Output projected points for debugging
  vtkSmartPointer<vtkVertexGlyphFilter> vertexGlyphFilter =
    vtkSmartPointer<vtkVertexGlyphFilter>::New();
  vertexGlyphFilter->AddInput(shiftTransformFilter->GetOutput());
  vertexGlyphFilter->Update();

  vtkSmartPointer<vtkXMLPolyDataWriter> writer =
    vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName("shiftedTransformedProjectedPoints.vtp");
  writer->SetInputConnection(vertexGlyphFilter->GetOutputPort());
  writer->Write();
  }

  // Create an image into which to project the points
  typedef itk::Image<itk::CovariantVector<unsigned char, 3>, 2> RGBImageType;
  RGBImageType::Pointer projectedImage = RGBImageType::New();
  projectedImage->SetRegions(this->FullImage->GetLargestPossibleRegion());
  projectedImage->Allocate();
  itk::CovariantVector<unsigned char, 3> black;
  black[0] = 0;
  black[1] = 0;
  black[2] = 0;
  projectedImage->FillBuffer(black);

  double newBounds[6];
  shiftTransformFilter->GetOutput()->GetBounds(newBounds);

  float binSize[2];
  binSize[0] = newBounds[1] / static_cast<float>(projectedImage->GetLargestPossibleRegion().GetSize()[0]); // xmax / width
  binSize[1] = newBounds[3] / static_cast<float>(projectedImage->GetLargestPossibleRegion().GetSize()[1]); // ymax / height

  for(unsigned int i = 0; i < static_cast<unsigned int>(shiftTransformFilter->GetOutput()->GetNumberOfPoints()); i++)
    {
    // Divide the 3D coordinates by the number of pixels (aka bins) to get the image coordinate
    double p[3];
    shiftTransformFilter->GetOutput()->GetPoint(i,p);
    itk::Index<2> index;
    index[0] = floor(p[0] / binSize[0]);
    index[1] = floor(p[1] / binSize[1]);
    //std::cout << "Index: " << index << std::endl;

    unsigned char pointColor[3];
    vtkUnsignedCharArray::SafeDownCast(shiftTransformFilter->GetOutput()->GetPointData()->GetScalars())->GetTupleValue(i, pointColor);

    itk::CovariantVector<unsigned char, 3> color;
    color[0] = pointColor[0];
    color[1] = pointColor[1];
    color[2] = pointColor[2];

    //std::cout << "color: " << color << std::endl;

    projectedImage->SetPixel(index, color);
    }

  typedef  itk::ImageFileWriter<RGBImageType> RGBWriterType;
  RGBWriterType::Pointer writer = RGBWriterType::New();
  writer->SetFileName(filename);
  writer->SetInput(projectedImage);
  writer->Update();
}
