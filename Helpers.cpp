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
#include <vtkMath.h>
#include <vtkPolyData.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkXMLPolyDataWriter.h>

// ITK
#include "itkImageRegionIterator.h"

namespace Helpers
{

void OutputPolyData(vtkPolyData* const points, const std::string& filename)
{
  // Output projected points for debugging
  vtkSmartPointer<vtkVertexGlyphFilter> vertexGlyphFilter =
    vtkSmartPointer<vtkVertexGlyphFilter>::New();
  vertexGlyphFilter->SetInputData(points);
  vertexGlyphFilter->Update();

  vtkSmartPointer<vtkXMLPolyDataWriter> writer =
    vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName(filename.c_str());
  writer->SetInputConnection(vertexGlyphFilter->GetOutputPort());
  writer->Write();
}

void ITKRGBImageToVTKImage(const itk::Image<itk::RGBPixel<unsigned char>, 2>* const image,
                           vtkImageData* const outputImage)
{
  typedef itk::Image<itk::RGBPixel<unsigned char>, 2> RGBImageType;
  // Setup and allocate the VTK image
  //outputImage->SetNumberOfScalarComponents(3);
  //outputImage->SetScalarTypeToUnsignedChar();
  outputImage->SetDimensions(image->GetLargestPossibleRegion().GetSize()[0],
                             image->GetLargestPossibleRegion().GetSize()[1],
                             1);

  //outputImage->AllocateScalars();
  outputImage->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

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

void ITKVectorImageToRGBImage(const itk::VectorImage<float, 2>* const image,
                              itk::Image<itk::RGBPixel<unsigned char>, 2>* const outputImage)
{
  if(image->GetNumberOfComponentsPerPixel() < 3)
  {
    throw std::runtime_error("Cannot convert a vector image with < 3 channels to an RGB image!");
  }

  outputImage->SetRegions(image->GetLargestPossibleRegion());
  outputImage->Allocate();

  typedef itk::Image<itk::RGBPixel<unsigned char>, 2> RGBImageType;
  typedef itk::VectorImage<float, 2> VectorImageType;

  itk::ImageRegionConstIteratorWithIndex<VectorImageType> imageIterator(image, image->GetLargestPossibleRegion());

  while(!imageIterator.IsAtEnd())
    {
    RGBImageType::PixelType outputPixel;
    outputPixel.SetRed(imageIterator.Get()[0]);
    outputPixel.SetGreen(imageIterator.Get()[1]);
    outputPixel.SetBlue(imageIterator.Get()[2]);
    outputImage->SetPixel(imageIterator.GetIndex(), outputPixel);

    ++imageIterator;
    }
}

void ITKRGBImageToVectorImage(const itk::Image<itk::RGBPixel<unsigned char>, 2>* const image,
                              itk::VectorImage<float, 2>* const outputImage)
{
  outputImage->SetRegions(image->GetLargestPossibleRegion());
  outputImage->SetNumberOfComponentsPerPixel(3);
  outputImage->Allocate();

  typedef itk::Image<itk::RGBPixel<unsigned char>, 2> RGBImageType;
  typedef itk::VectorImage<float, 2> VectorImageType;
  
  itk::ImageRegionConstIteratorWithIndex<RGBImageType> imageIterator(image, image->GetLargestPossibleRegion());

  while(!imageIterator.IsAtEnd())
    {
    // Can't get a reference of a pixel from a VectorImage apparently?
    //VectorImageType::PixelType& outputPixel = outputImage->GetPixel(imageIterator.GetIndex());
    VectorImageType::PixelType outputPixel = outputImage->GetPixel(imageIterator.GetIndex());
    outputPixel[0] = imageIterator.Get().GetRed();
    outputPixel[1] = imageIterator.Get().GetGreen();
    outputPixel[2] = imageIterator.Get().GetBlue();
    outputImage->SetPixel(imageIterator.GetIndex(), outputPixel);
  
    ++imageIterator;
    }
}

unsigned int NumberOfUniquePoints(vtkPoints* const points, const float tolerance)
{
  if(points->GetNumberOfPoints() <= 1)
  {
    return points->GetNumberOfPoints();
  }
  
  double currentPoint[3];
  points->GetPoint(0, currentPoint);
  unsigned int numberOfUniquePoints = 1;
  
  double p[3];

  for(vtkIdType i = 1; i < points->GetNumberOfPoints(); ++i)
    {
    points->GetPoint(i,p);
    double distance = sqrt(vtkMath::Distance2BetweenPoints(currentPoint, p));
    if(distance > tolerance)
      {
      points->GetPoint(i,currentPoint);
      numberOfUniquePoints++;
      }
    }
  return numberOfUniquePoints;
}


/******************************//**
 * \brief Conversion from cartesian to spherical coordinates (if the
 *        resulting radius is 0 also the azimutal and inclination
 *        angles get set to 0).
 *
 * @param r      Reference to return the radius
 * @param theta  Reference to return the inclination angle (in radian)
 * @param phi    Reference to return the azimutal angle (in radian)
 * @param x      x coordinate
 * @param y      y coordinate
 * @param z      z coordinate
 ******************************/
void
cartesianToSpherical( double & r,
                      double & theta,
                      double & phi,
                      double   x,
                      double   y,
                      double   z )
{
    if ( ( r = sqrt( x * x + y * y + z * z ) ) != 0.0 )
    {
        theta = acos( z / r );
        phi   = atan2( y, x );
    }
    else
        theta = phi = 0.0;
}

/******************************//**
 * \brief Conversion from spherical to cartesian coordinates
 *
 * @param x      Reference to return the x coordinate
 * @param y      Reference to return the y coordinate
 * @param z      Reference to return the z coordinate
 * @param r      Radius (must be non-negative)
 * @param theta  Inclination angle (in radian)
 * @param phi    Azimutal angle (in radian)
 ******************************/
void
sphericalToCartesian( double & x,
                      double & y,
                      double & z,
                      double   r,
                      double   theta,
                      double   phi )
{
        if ( r < 0.0 )
                throw "Negative radius in sphericalToCartesian()";
    x = r * sin( theta ) * cos( phi );
    y = r * sin( theta ) * sin( phi );
    z = r * cos( theta );
} 

void PrintSpherical(double x, double y, double z)
{
  double r, theta, phi;
  cartesianToSpherical(r, theta, phi, x, y, z );
  std::cout << "r: " << r << " theta: " << theta << " phi: " << phi << std::endl;
}

}; // end Helpers namespace