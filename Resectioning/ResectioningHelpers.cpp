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

#include "ResectioningHelpers.h"

// ITK
#include "itkImageRegionIterator.h"
#include "itkVectorMagnitudeImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"

// VTK
#include <vtkCell.h>
#include <vtkIdList.h>
#include <vtkKdTree.h>
#include <vtkMath.h>
#include <vtkPolyData.h>
#include <vtkStructuredGrid.h>

namespace ResectioningHelpers
{

// Convert a vector ITK image to a VTK image for display
void ITKImagetoVTKImage(const FloatVectorImageType* const image, vtkImageData* const outputImage)
{
  //std::cout << "ITKImagetoVTKImage()" << std::endl;
  if(image->GetNumberOfComponentsPerPixel() >= 3)
    {
    ITKImagetoVTKRGBImage(image, outputImage);
    }
  else
    {
    ITKImagetoVTKMagnitudeImage(image, outputImage);
    }
}

// Convert a vector ITK image to a VTK image for display
void ITKImagetoVTKRGBImage(const FloatVectorImageType* const image, vtkImageData* const outputImage)
{
  // This function assumes an ND (with N>3) image has the first 3 channels as
  // RGB and extra information in the remaining channels.
  
  std::cout << "ITKImagetoVTKRGBImage()" << std::endl;
  if(image->GetNumberOfComponentsPerPixel() < 3)
    {
    std::stringstream ss;
    ss << "The input image has " << image->GetNumberOfComponentsPerPixel()
       << " components, but at least 3 are required.";
    throw std::runtime_error(ss.str());
    }

  // Setup and allocate the image data
  //outputImage->SetNumberOfScalarComponents(3);
  //outputImage->SetScalarTypeToUnsignedChar();
  outputImage->SetDimensions(image->GetLargestPossibleRegion().GetSize()[0],
                             image->GetLargestPossibleRegion().GetSize()[1],
                             1);

  //outputImage->AllocateScalars();
  outputImage->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

  // Copy all of the input image pixels to the output image
  itk::ImageRegionConstIteratorWithIndex<FloatVectorImageType> imageIterator(image,
                                                                             image->GetLargestPossibleRegion());
  imageIterator.GoToBegin();

  while(!imageIterator.IsAtEnd())
    {
    unsigned char* pixel = static_cast<unsigned char*>(outputImage->GetScalarPointer(imageIterator.GetIndex()[0],
                                                                                     imageIterator.GetIndex()[1],0));
    for(unsigned int component = 0; component < 3; ++component)
      {
      unsigned char value = static_cast<unsigned char>(imageIterator.Get()[component]);
      //std::cout << "Value: " << static_cast<int>(value) << std::endl;
      pixel[component] = value;
      }

    ++imageIterator;
    }
}


// Convert a vector ITK image to a VTK image for display
void ITKImagetoVTKMagnitudeImage(const FloatVectorImageType* const image, vtkImageData* const outputImage)
{
  //std::cout << "ITKImagetoVTKMagnitudeImage()" << std::endl;
  // Compute the magnitude of the ITK image
  typedef itk::VectorMagnitudeImageFilter<
                  FloatVectorImageType, FloatScalarImageType >  VectorMagnitudeFilterType;

  // Create and setup a magnitude filter
  VectorMagnitudeFilterType::Pointer magnitudeFilter = VectorMagnitudeFilterType::New();
  magnitudeFilter->SetInput( image );
  magnitudeFilter->Update();

  // Rescale and cast for display
  typedef itk::RescaleIntensityImageFilter<
                  FloatScalarImageType, UnsignedCharScalarImageType > RescaleFilterType;

  RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
  rescaleFilter->SetOutputMinimum(0);
  rescaleFilter->SetOutputMaximum(255);
  rescaleFilter->SetInput( magnitudeFilter->GetOutput() );
  rescaleFilter->Update();

  // Setup and allocate the VTK image
  //outputImage->SetNumberOfScalarComponents(1);
  //outputImage->SetScalarTypeToUnsignedChar();
  outputImage->SetDimensions(image->GetLargestPossibleRegion().GetSize()[0],
                             image->GetLargestPossibleRegion().GetSize()[1],
                             1);

  //outputImage->AllocateScalars();
  outputImage->AllocateScalars(VTK_UNSIGNED_CHAR, 1);

  // Copy all of the scaled magnitudes to the output image
  itk::ImageRegionConstIteratorWithIndex<UnsignedCharScalarImageType> imageIterator(rescaleFilter->GetOutput(), rescaleFilter->GetOutput()->GetLargestPossibleRegion());
  imageIterator.GoToBegin();

  while(!imageIterator.IsAtEnd())
    {
    unsigned char* pixel = static_cast<unsigned char*>(outputImage->GetScalarPointer(imageIterator.GetIndex()[0],
                                                                                     imageIterator.GetIndex()[1],0));
    pixel[0] = imageIterator.Get();

    ++imageIterator;
    }
}

float ComputeAverageSpacing(vtkPoints* const points, unsigned int numberOfPointsToUse)
{
  // Compute the average spacing between pairs of closest points in the data.
  // If 'numberOfPointsToUse' is anything positive, only the first 'numberOfPointsToUse'
  // in the data set will be used in the computation. This is an option because often
  // in very large datasets (~1M points) a small sample is enough to determine the average
  // spacing relatively accurately.
  
  if(numberOfPointsToUse == 0)
    {
    numberOfPointsToUse = points->GetNumberOfPoints();
    }
  if(numberOfPointsToUse > static_cast<unsigned int>(points->GetNumberOfPoints()))
    {
    numberOfPointsToUse = points->GetNumberOfPoints();
    }

  float sumOfDistances = 0.;
  //Create the tree
  vtkSmartPointer<vtkKdTree> pointTree = vtkSmartPointer<vtkKdTree>::New();
  pointTree->BuildLocatorFromPoints(points);

  for(vtkIdType i = 0; i < static_cast<vtkIdType>(numberOfPointsToUse); ++i)
    {
    // Get the coordinates of the current point
    double queryPoint[3];
    points->GetPoint(i,queryPoint);
  
    // Find the 2 closest points (the first closest will be exactly the query point)
    vtkSmartPointer<vtkIdList> result = vtkSmartPointer<vtkIdList>::New();
    pointTree->FindClosestNPoints(2, queryPoint, result);
  
    double closestPoint[3];
    points->GetPoint(result->GetId(1), closestPoint);
      
    float squaredDistance = vtkMath::Distance2BetweenPoints(queryPoint, closestPoint);

    // Take the square root to get the Euclidean distance between the points.
    float distance = sqrt(squaredDistance);
    
    sumOfDistances += distance;
    }
    
  float averageDistance = sumOfDistances / static_cast<float>(points->GetNumberOfPoints());
  
  return averageDistance;
}


void StructuredGridToPolyData(vtkStructuredGrid* const structuredGrid, vtkPolyData* const polyData)
{
  polyData->SetPoints(structuredGrid->GetPoints());

  polyData->Allocate();

  for(vtkIdType cellId = 0; cellId < structuredGrid->GetNumberOfCells(); ++cellId)
    {
    vtkCell* cell = structuredGrid->GetCell(cellId);
    //polyData->InsertNextCell(cell->GetCellType(), cell->GetNumberOfPoints(), cell->GetPoints());
    vtkIdType pointIds[cell->GetNumberOfPoints()];
    for(vtkIdType pointId = 0; pointId < cell->GetNumberOfPoints(); ++pointId)
      {
      pointIds[pointId] = cell->GetPointId(pointId);
      }
    polyData->InsertNextCell(cell->GetCellType(), cell->GetNumberOfPoints(), pointIds);

    //std::cout << "cell " << cellId << " has " << cell->GetNumberOfPoints() << " points." << std::endl;
    }
}

} // end namespace
