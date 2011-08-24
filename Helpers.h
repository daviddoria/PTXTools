#ifndef HELPERS_H
#define HELPERS_H

// VTK
#include "vtkSmartPointer.h"

// ITK
#include "itkImageRegionIterator.h"

class vtkPolyData;

namespace Helpers
{
  void OutputPolyData(vtkSmartPointer<vtkPolyData> points, std::string filename);

  template<typename TImage>
  void SetAllPixels(typename TImage::Pointer image, typename TImage::PixelType pixel);

  template<typename TImage>
  void DeepCopy(typename TImage::Pointer input, typename TImage::Pointer output);
};

#include "Helpers.txx"

#endif
