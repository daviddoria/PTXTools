/*=========================================================================
 *
 *  Copyright David Doria 2011 daviddoria@gmail.com
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

#ifndef ResectioningHelpers_H
#define ResectioningHelpers_H

// ITK
#include "itkImage.h"
#include "itkIndex.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"

// VTK
#include <vtkImageData.h>
#include <vtkPoints.h>
#include <vtkProp.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
class vtkStructuredGrid;
class vtkPolyData;
// Custom
#include "Types.h"

namespace ResectioningHelpers
{
  
template<typename T>
void RemoveAllActors(std::vector<T> actors, vtkRenderer* const renderer);
//void RemoveAllActors(const std::vector<T>& actors, vtkRenderer* const renderer);

// This function simply drives ITKImagetoVTKRGBImage or ITKImagetoVTKMagnitudeImage
void ITKImagetoVTKImage(const FloatVectorImageType* const image, vtkImageData* const outputImage); 
void ITKImagetoVTKRGBImage(const FloatVectorImageType* const image, vtkImageData* const outputImage);
void ITKImagetoVTKMagnitudeImage(const FloatVectorImageType* const image, vtkImageData* const outputImage);

/** Compute the average distance between neighboring points. 'numberOfPointsToUse'
 * is not const because if it is 0, then it is set to the full number of points in the data set. */
float ComputeAverageSpacing(vtkPoints* const points, unsigned int numberOfPointsToUse);

void StructuredGridToPolyData(vtkStructuredGrid* const structuredGrid, vtkPolyData* const polyData);

}


#include "ResectioningHelpers.hpp"

#endif
