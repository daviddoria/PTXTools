/*
Copyright (C) 2010 David Doria, daviddoria@gmail.com

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

#include "Helpers.h"

// STL
#include <string>

// VTK
#include <vtkPolyData.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkXMLPolyDataWriter.h>

// ITK
#include "itkImageRegionIterator.h"

namespace Helpers
{

void OutputPolyData(vtkSmartPointer<vtkPolyData> points, std::string filename)
{
  // Output projected points for debugging
  vtkSmartPointer<vtkVertexGlyphFilter> vertexGlyphFilter =
    vtkSmartPointer<vtkVertexGlyphFilter>::New();
  vertexGlyphFilter->AddInput(points);
  vertexGlyphFilter->Update();

  vtkSmartPointer<vtkXMLPolyDataWriter> writer =
    vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName(filename.c_str());
  writer->SetInputConnection(vertexGlyphFilter->GetOutputPort());
  writer->Write();
}

void ITKRGBImageToVTKImage(const itk::Image<itk::RGBPixel<unsigned char>, 2>::Pointer image, vtkImageData* outputImage)
{
  typedef itk::Image<itk::RGBPixel<unsigned char>, 2> RGBImageType;
  // Setup and allocate the VTK image
  outputImage->SetNumberOfScalarComponents(3);
  outputImage->SetScalarTypeToUnsignedChar();
  outputImage->SetDimensions(image->GetLargestPossibleRegion().GetSize()[0],
                             image->GetLargestPossibleRegion().GetSize()[1],
                             1);

  outputImage->AllocateScalars();

  // Copy all of the scaled magnitudes to the output image
  itk::ImageRegionConstIteratorWithIndex<RGBImageType> imageIterator(image, image->GetLargestPossibleRegion());
  imageIterator.GoToBegin();

  while(!imageIterator.IsAtEnd())
    {
    unsigned char* pixel = static_cast<unsigned char*>(outputImage->GetScalarPointer(imageIterator.GetIndex()[0],
                                                                                     imageIterator.GetIndex()[1],0));
    pixel[0] = imageIterator.Get().GetRed();
    pixel[1] = imageIterator.Get().GetGreen();
    pixel[2] = imageIterator.Get().GetBlue();

    ++imageIterator;
    }

  outputImage->Modified();
}

}; // end Helpers namespace