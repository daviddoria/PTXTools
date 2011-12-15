#include "PTXImage.h"

#include "itkImageFileReader.h"
#include "itkImageRegionConstIterator.h"

#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkXMLPolyDataWriter.h"
#include <vtkVertexGlyphFilter.h>

#include <string>

// This program takes a binary image. It extracts the points from the PTX image which correspond to non-zero pixels in the binary image.

int main (int argc, char *argv[])
{
  // Verify arguments
  if(argc != 4)
    {
    std::cerr << "Required arguments: InputFilename(ptx) MaskBoundaryFileName(png) OutputFilename(vtp)" << std::endl;
    return EXIT_FAILURE;
    }

  // Parse arguments
  std::string ptxFileName = argv[1];
  std::string maskFileName = argv[2];
  std::string outputPointsFileName = argv[3];

  std::cout << "PTX file: " << ptxFileName << std::endl;
  std::cout << "Mask file: " << maskFileName << std::endl;
  std::cout << "Output points file: " << outputPointsFileName << std::endl;

  // Read the mask image
  typedef itk::Image< unsigned char, 2 >         MaskImageType;
  typedef itk::ImageFileReader<MaskImageType>    MaskReaderType;

  MaskReaderType::Pointer maskReader = MaskReaderType::New();
  maskReader->SetFileName(maskFileName);
  maskReader->Update();

  PTXImage ptxImage;
  ptxImage.ReadFile(ptxFileName);

  vtkSmartPointer<vtkPoints> points =
    vtkSmartPointer<vtkPoints>::New();

  itk::ImageRegionConstIterator<MaskImageType> imageIterator(maskReader->GetOutput(), maskReader->GetOutput()->GetLargestPossibleRegion());

  while(!imageIterator.IsAtEnd())
    {
    if(imageIterator.Get()) // the pixel is non-zero
      {
      // Get the point from the ptx image
      PTXPixel pixel = ptxImage.GetFullImage()->GetPixel(imageIterator.GetIndex());

      // Store it in the list of points
      points->InsertNextPoint ( pixel.X, pixel.Y, pixel.Z );
      }
    ++imageIterator;
    }

  vtkSmartPointer<vtkPolyData> polydata =
    vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints(points);

  vtkSmartPointer<vtkVertexGlyphFilter> vertexGlyphFilter =
    vtkSmartPointer<vtkVertexGlyphFilter>::New();
  vertexGlyphFilter->AddInput(polydata);
  vertexGlyphFilter->Update();

  // Write the file
  vtkSmartPointer<vtkXMLPolyDataWriter> writer =
    vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName(outputPointsFileName.c_str());
  writer->SetInputConnection(vertexGlyphFilter->GetOutputPort());
  writer->Write();

  return EXIT_SUCCESS;
}
