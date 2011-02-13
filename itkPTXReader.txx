#ifndef __itkPTXReader_txx
#define __itkPTXReader_txx

#include "itkPTXReader.h"

#include "itkObjectFactory.h"
#include "itkImageRegionIterator.h"
#include "itkImageFileWriter.h"

#include <fstream>
#include <sstream>
#include <string>

namespace itk
{
typedef itk::Image<itk::CovariantVector< float, 8>, 2> FullImageType;

FullImageType::Pointer PTXReader::GetFullImage()
{
  return this->FullImage;
}

void PTXReader::WriteRGBImage(std::string filename)
{

  typedef itk::Image<itk::CovariantVector<unsigned char, 3>, 2> RGBImageType;
  RGBImageType::Pointer rgbImage = RGBImageType::New();
  rgbImage->SetRegions(this->FullImage->GetLargestPossibleRegion());
  rgbImage->Allocate();

  itk::ImageRegionIterator<RGBImageType> rgbImageIterator(rgbImage,rgbImage->GetLargestPossibleRegion());
  rgbImageIterator.GoToBegin();

  itk::ImageRegionIterator<FullImageType> fullImageIterator(this->FullImage, this->FullImage->GetLargestPossibleRegion());
  fullImageIterator.GoToBegin();

  while(!rgbImageIterator.IsAtEnd())
    {
    itk::CovariantVector<float, 8> fullPixel = fullImageIterator.Get();

    itk::CovariantVector<unsigned char, 3> rgbPixel;
    rgbPixel[0] = round(fullPixel[4]);
    rgbPixel[1] = round(fullPixel[5]);
    rgbPixel[2] = round(fullPixel[6]);

    rgbImageIterator.Set(rgbPixel);

    ++rgbImageIterator;
    ++fullImageIterator;
    }

  typedef  itk::ImageFileWriter< RGBImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(filename);
  writer->SetInput(rgbImage);
  writer->Update();

}


void PTXReader::GenerateData()
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

  itk::Image<itk::CovariantVector<float, 8>, 2 >::Pointer output = this->GetOutput();

  itk::Index<2> start;
  start[0] = 0;
  start[1] = 0;

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

      itk::CovariantVector<float, 8> pixel;
      pixel[0] = coordinate[0];
      pixel[1] = coordinate[1];
      pixel[2] = coordinate[2];
      pixel[3] = intensity;
      pixel[4] = color[0];
      pixel[5] = color[1];
      pixel[6] = color[2];
      pixel[7] = point.GetNorm();

      output->SetPixel(pixelIndex, pixel);
      }//end phi for loop
    }// end theta for loop

  // Close the input file
  infile.close();

  this->FullImage = FullImageType::New();
  this->FullImage->Graft(output);

}

}// end namespace


#endif
