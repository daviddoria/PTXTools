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

  /** Write a PolyData to a file */
  void OutputPolyData(vtkPolyData* const points, const std::string& filename);

  /** Convert an ITK image to a VTK image */
  void ITKRGBImageToVTKImage(const itk::Image<itk::RGBPixel<unsigned char>, 2>* const image, vtkImageData* const outputImage);

  ///// Templated functions /////
  /** Set every pixel in an image to 'pixel' */
  template<typename TImage>
  void SetAllPixelsToValue(TImage* const image, typename TImage::PixelType& pixel);

  /** Deep copy an image. */
  template<typename TImage>
  void DeepCopy(const TImage* const input, TImage* const output);

  /** Rescale and convert an ITK image to a VTK image*/
  template <typename TImage>
  void ITKScalarImageToScaledVTKImage(const TImage* const image, vtkImageData* const outputImage);

  /** Write an image to a file. */
  template<typename TImage>
  void WriteImage(const TImage* const image, const std::string& filename);

  /** Get the nearest valid pixel location. */
  template<typename TImage>
  itk::Index<2> GetNearestValidPixel(const TImage* const image, const itk::Index<2>& queryPixel);

};

#include "Helpers.hpp"

#endif
