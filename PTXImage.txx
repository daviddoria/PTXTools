#include "PTXImage.h"

// ITK
#include "itkCovariantVector.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"

// VTK
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkUnsignedCharArray.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkXMLPolyDataWriter.h>

// STL
#include <fstream>
#include <sstream>
#include <string>

PTXImage::PTXImage()
{
  this->FullImage = FullImageType::New();
}

void PTXImage::ReadFile(std::string filename)
{
  std::ifstream infile;
  infile.open(filename.c_str());
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

      PTXPixel pixel;
      pixel.X = coordinate[0];
      pixel.Y = coordinate[1];
      pixel.Z = coordinate[2];
      pixel.Intensity = intensity;
      pixel.R = color[0];
      pixel.G = color[1];
      pixel.B = color[2];
      //std::cout << pixel << std::endl;

      this->FullImage->SetPixel(pixelIndex, pixel);
      }//end phi for loop
    }// end theta for loop

  // Close the input file
  infile.close();

}

void PTXImage::WriteRGBImage(std::string filename)
{
  typedef itk::Image<itk::CovariantVector<unsigned char, 3>, 2> RGBImageType;

  RGBImageType::Pointer rgbImage = RGBImageType::New();
  rgbImage->SetRegions(this->FullImage->GetLargestPossibleRegion());
  rgbImage->Allocate();

  itk::ImageRegionIterator<RGBImageType> rgbImageIterator(rgbImage, rgbImage->GetLargestPossibleRegion());
  rgbImageIterator.GoToBegin();

  itk::ImageRegionIterator<FullImageType> fullImageIterator(this->FullImage, this->FullImage->GetLargestPossibleRegion());
  fullImageIterator.GoToBegin();

  while(!rgbImageIterator.IsAtEnd())
    {
    PTXPixel fullPixel = fullImageIterator.Get();

    itk::CovariantVector<unsigned char, 3> rgbPixel;
    rgbPixel[0] = fullPixel.R;
    rgbPixel[1] = fullPixel.G;
    rgbPixel[2] = fullPixel.B;

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


void PTXImage::WritePointCloud(std::string filename)
{
  vtkSmartPointer<vtkPoints> points =
    vtkSmartPointer<vtkPoints>::New();

  vtkSmartPointer<vtkUnsignedCharArray> colors =
    vtkSmartPointer<vtkUnsignedCharArray>::New();
  colors->SetNumberOfComponents(3);
  colors->SetName("Colors");

  itk::ImageRegionConstIteratorWithIndex<FullImageType> imageIterator(this->FullImage, this->FullImage->GetLargestPossibleRegion());

  while(!imageIterator.IsAtEnd())
    {
    PTXPixel pixel = imageIterator.Get();
    unsigned char rgb[3];
    rgb[0] = pixel.R;
    rgb[1] = pixel.G;
    rgb[2] = pixel.B;
    colors->InsertNextTupleValue(rgb);

    points->InsertNextPoint(pixel.X, pixel.Y, pixel.Z);

    ++imageIterator;
    }

  vtkSmartPointer<vtkPolyData> polydata =
    vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints(points);
  polydata->GetPointData()->SetScalars(colors);

  vtkSmartPointer<vtkVertexGlyphFilter> vertexGlyphFilter =
    vtkSmartPointer<vtkVertexGlyphFilter>::New();
  vertexGlyphFilter->AddInput(polydata);
  vertexGlyphFilter->Update();

  vtkSmartPointer<vtkXMLPolyDataWriter> writer =
    vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName(filename.c_str());
  writer->SetInputConnection(vertexGlyphFilter->GetOutputPort());
  writer->Write();
}