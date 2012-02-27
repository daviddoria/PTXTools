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

#ifndef HELPERS_H
#define HELPERS_H

// ITK
#include "itkImage.h"
#include "itkRGBPixel.h"

// VTK
#include "vtkSmartPointer.h"
class vtkImageData;
class vtkPolyData;

namespace Helpers
{
  ///// Non-templated functions /////
  void OutputPolyData(vtkSmartPointer<vtkPolyData> points, std::string filename);
  void ITKRGBImageToVTKImage(const itk::Image<itk::RGBPixel<unsigned char>, 2>::Pointer image, vtkImageData* outputImage);

  ///// Templated functions /////
  template<typename TImage>
  void SetAllPixels(typename TImage::Pointer image, typename TImage::PixelType pixel);

  template<typename TImage>
  void DeepCopy(const TImage* const input, TImage* const output);

  template <typename TImage>
  void ITKScalarImageToScaledVTKImage(const typename TImage::Pointer image, vtkImageData* outputImage);

  template<typename TImage>
  void WriteImage(const TImage* const image, const std::string& filename);

};

#include "Helpers.hxx"

#endif
