// STL
#include <iostream>

// VTK
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataReader.h>

// ITK
#include "itkImage.h"
#include "itkImageFileWriter.h"

// Custom
#include "Helpers.h"
#include "PTXImage.h"

int main(int argc, char *argv[])
{
  if(argc < 3)
    {
    std::cerr << "Required arguments: pointcloud.vtp output.mha" << std::endl;
    return EXIT_FAILURE;
    }
  // Parse arguments
  std::string inputFileName = argv[1];
  std::string outputFileName = argv[2];

  // Output arguments
  std::cout << "inputFileName : " << inputFileName << std::endl;
  std::cout << "outputFileName : " << outputFileName << std::endl;

  // Read the VTP file
  vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
  reader->SetFileName(inputFileName.c_str());
  reader->Update();

  typedef itk::VectorImage<float, 2> FloatVectorImageType;
  FloatVectorImageType::Pointer normalsImage = FloatVectorImageType::New();
  normalsImage->SetNumberOfComponentsPerPixel(3);
  
  vtkIntArray* imageSizeArray = vtkIntArray::SafeDownCast ( reader->GetOutput()->GetFieldData()->GetArray ( "ImageSize" ) );
  int imageSize[2];
  imageSizeArray->GetTupleValue(0, imageSize);
  itk::Index<2> corner;
  corner.Fill(0);
  
  itk::Size<2> size;
  size[0] = imageSize[0];
  size[1] = imageSize[1];
  
  itk::ImageRegion<2> region(corner,size);
  normalsImage->SetRegions(region);
  normalsImage->Allocate();
  
  vtkIntArray* originalPixelArray = vtkIntArray::SafeDownCast ( reader->GetOutput()->GetPointData()->GetArray ( "OriginalPixel" ) );
  
  if(!originalPixelArray)
    {
    std::cerr << "VTP file does not contain OriginalPixel array!" << std::endl;
    return EXIT_FAILURE;
    }

  vtkFloatArray* normalsArray = vtkFloatArray::SafeDownCast ( reader->GetOutput()->GetPointData()->GetNormals() );
  
  if(!normalsArray)
    {
    std::cerr << "VTP file does not contain point normals!" << std::endl;
    return EXIT_FAILURE;
    }

  for(vtkIdType pointId = 0; pointId < reader->GetOutput()->GetNumberOfPoints(); pointId++)
    {
    int originalPixel[2];
    originalPixelArray->GetTupleValue(pointId, originalPixel);
    itk::Index<2> index;
    index[0] = originalPixel[0];
    index[1] = originalPixel[1];
  
    FloatVectorImageType::PixelType pixel;
    pixel.SetSize(3);
  
    float normal[3];
    normalsArray->GetTupleValue(pointId, normal);
    
    pixel[0] = normal[0];
    pixel[1] = normal[1];
    pixel[2] = normal[2];
    
    normalsImage->SetPixel(index, pixel);
    }
  
  typedef itk::ImageFileWriter<FloatVectorImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(outputFileName);
  writer->SetInput(normalsImage);
  writer->Update();

  return EXIT_SUCCESS;
}
