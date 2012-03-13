/*
Copyright (C) 2011 David Doria, daviddoria@gmail.com

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

#include "itkImage.h"

template <typename TImage>
ImageLayer<TImage>::ImageLayer()
{
  this->Image = TImage::New();
  itk::Index<2> corner;
  corner.Fill(0);
  itk::Size<2> size;
  size.Fill(0);
  itk::ImageRegion<2> region(corner,size);
  this->Image->SetRegions(region);
  this->Image->Allocate();
  
  this->ImageData = vtkSmartPointer<vtkImageData>::New();
  this->ImageData->SetDimensions(0,0,0);
  //this->ImageData->AllocateScalars();

  this->ImageSlice = vtkSmartPointer<vtkImageSlice>::New();
  this->ImageSliceMapper = vtkSmartPointer<vtkImageSliceMapper>::New();
  this->ImageSliceMapper->SetInputData(this->ImageData);
  this->ImageSlice->SetMapper(this->ImageSliceMapper);
  this->ImageSlice->VisibilityOff();

}