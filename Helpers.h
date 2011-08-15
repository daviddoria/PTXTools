#ifndef HELPERS_H
#define HELPERS_H

#include "vtkSmartPointer.h"

#include "itkImageRegionIterator.h"

class vtkPolyData;

namespace Helpers
{
  void OutputPolyData(vtkSmartPointer<vtkPolyData> points, std::string filename);

  /** Copy the input to the output*/
  template<typename TImage>
  void DeepCopy(typename TImage::Pointer input, typename TImage::Pointer output)
  {
    output->SetRegions(input->GetLargestPossibleRegion());
    output->Allocate();

    itk::ImageRegionConstIterator<TImage> inputIterator(input, input->GetLargestPossibleRegion());
    itk::ImageRegionIterator<TImage> outputIterator(output, output->GetLargestPossibleRegion());

    while(!inputIterator.IsAtEnd())
      {
      outputIterator.Set(inputIterator.Get());
      ++inputIterator;
      ++outputIterator;
      }
  }
};

#endif