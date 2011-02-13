#ifndef __itkPTXRGBDIReader_txx
#define __itkPTXRGBDIReader_txx

#include "itkPTXRGBDIReader.h"

#include "itkObjectFactory.h"
#include "itkImageRegionIterator.h"

#include <fstream>
#include <sstream>

namespace itk
{

void PTXRGBDIReader::GenerateData()
{
  std::ifstream infile;
  infile.open(this->m_FileName.c_str());
  if(!infile)
    {
    std::cout << "Could not open ptx file " << this->m_FileName << "!" << std::endl;
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

  // Skip 8 lines (identity matrices)
  for(int i = 0; i < 8; i++)
    {
    getline(infile, line);
    }

  std::cout << "PhiPoints: " << numberOfPhiPoints << std::endl;
  std::cout << "ThetaPoints: " << numberOfThetaPoints << std::endl;

  typedef itk::Image<itk::CovariantVector<float, 5>, 2 > ImageType;
  ImageType::Pointer output = this->GetOutput();

  itk::Index<2> start;
  start.Fill(0);

  itk::Size<2> size;
  size[0] = numberOfThetaPoints;
  size[1] = numberOfPhiPoints;

  itk::ImageRegion<2> region;
  region.SetSize(size);
  region.SetIndex(start);

  output->SetRegions(region);
  output->Allocate();

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
      float color[3];

      std::stringstream parsedLine(line);
      parsedLine >> coordinate[0] >> coordinate[1] >> coordinate[2] >> intensity
                 >> color[0] >> color[1] >> color[2];

      itk::CovariantVector<float, 3> point;
      point[0] = coordinate[0];
      point[1] = coordinate[1];
      point[2] = coordinate[2];

      itk::CovariantVector<float, 5> pixel;
      pixel[0] = color[0];
      pixel[1] = color[1];
      pixel[2] = color[2];
      pixel[3] = point.GetNorm(); // Depth
      pixel[4] = intensity;

      output->SetPixel(pixelIndex, pixel);
      }//end phi for loop
    }// end theta for loop

  // Normalize components 3 and 4 (depth and intensity) to be between 0 and 255 (like the RGB components) so that all components are in the same magnitude range

  // First, find the extremes
  itk::ImageRegionIterator<ImageType> imageIterator(output,region);

  float minDepth = std::numeric_limits<float>::max();
  float maxDepth = std::numeric_limits<float>::min();

  float minIntensity = std::numeric_limits<float>::max();
  float maxIntensity = std::numeric_limits<float>::min();

  imageIterator.GoToBegin();
  while(!imageIterator.IsAtEnd())
    {
    itk::CovariantVector<float, 5> val = imageIterator.Get();
    if(val[3] < minDepth)
      {
      minDepth = val[3];
      }
    if(val[3] > maxDepth)
      {
      maxDepth = val[3];
      }

    if(val[4] < minIntensity)
      {
      minIntensity = val[4];
      }
    if(val[4] > maxIntensity)
      {
      maxIntensity = val[4];
      }

    ++imageIterator;
    }

  std::cout << "minDepth: " << minDepth << std::endl;
  std::cout << "maxDepth: " << maxDepth << std::endl;

  std::cout << "minIntensity: " << minIntensity << std::endl;
  std::cout << "maxIntensity: " << maxIntensity << std::endl;
  // Perform the normalization
  imageIterator.GoToBegin();
  while(!imageIterator.IsAtEnd())
    {
    itk::CovariantVector<float, 5> val = imageIterator.Get();
    val[3] = 255. * (val[3] - minDepth)/(maxDepth - minDepth);
    val[4] = 255. * (val[4] - minIntensity)/(maxIntensity - minIntensity);

    imageIterator.Set(val);
    ++imageIterator;
    }

  // Close the input file
  infile.close();

}

}// end namespace


#endif
