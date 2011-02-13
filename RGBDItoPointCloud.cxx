#include "itkPTXReader.h"

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageRegionConstIteratorWithIndex.h>

#include <vtkXMLPolyDataWriter.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkUnsignedCharArray.h>
#include <vtkVertexGlyphFilter.h>

typedef itk::Image<itk::CovariantVector<float, 5> > RGBDIImageType;

void RGBDItoPointCloud(RGBDIImageType::Pointer rgbdiImage, vtkSmartPointer<vtkPolyData> pointCloud, float thetaStart, float phiStart, float thetaStep, float phiStep);

typedef itk::Image<itk::CovariantVector< float, 8>, 2> FullImageType;

int main (int argc, char *argv[])
{
  if(argc != 4)
    {
    std::cout << "Required arguments: InputFilename(mhd (RGBDI)) OriginalPTX OutputFilename(vtp)" << std::endl;
    return EXIT_FAILURE;
    }

  std::string inputFilename = argv[1];
  std::string originalPTX = argv[2];
  std::string outputFilename = argv[3];

  itk::ImageFileReader<RGBDIImageType>::Pointer reader = itk::ImageFileReader<RGBDIImageType>::New();
  reader->SetFileName(inputFilename);
  reader->Update();

  itk::PTXReader::Pointer ptxReader = itk::PTXReader::New();
  ptxReader->SetFileName(originalPTX);
  ptxReader->Update();


  itk::Index<2> index;
  index.Fill(0);
  FullImageType::Pointer fullImage = ptxReader->GetFullImage();
  FullImageType::PixelType corner = fullImage->GetPixel(index);

  float theta0 = atan(corner[0]/corner[1]);
  float phi0 = atan(corner[2]/sqrt(corner[0]*corner[0] + corner[1]*corner[1]));

  std::cout << "theta start: " << theta0 << std::endl;
  std::cout << "phi start: " << phi0 << std::endl;

  itk::Index<2> index1;
  index1[0] = 1;
  index1[1] = 0;
  FullImageType::PixelType pixel1 = fullImage->GetPixel(index1);
  float theta1 = atan(pixel1[0]/pixel1[1]);

  itk::Index<2> index2;
  index2[0] = 0;
  index2[1] = 1;
  FullImageType::PixelType pixel2 = fullImage->GetPixel(index2);
  float phi1 = atan(pixel2[2]/sqrt(pixel2[0]*pixel2[0] + pixel2[1]*pixel2[1]));

  float thetaStep = fabs(theta0 - theta1);
  float phiStep = fabs(phi0 - phi1);

  std::cout << "theta step: " << thetaStep << std::endl;
  std::cout << "phi step: " << phiStep << std::endl;

  vtkSmartPointer<vtkPolyData> pointCloud =
    vtkSmartPointer<vtkPolyData>::New();
  RGBDItoPointCloud(reader->GetOutput(), pointCloud, theta0, phi0, thetaStep, phiStep);

  vtkSmartPointer<vtkXMLPolyDataWriter> writer =
    vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName(outputFilename.c_str());
  writer->SetInputConnection(pointCloud->GetProducerPort());
  writer->Write();

  return EXIT_SUCCESS;
}

void RGBDItoPointCloud(RGBDIImageType::Pointer rgbdiImage, vtkSmartPointer<vtkPolyData> pointCloud, float thetaStart, float phiStart, float thetaStep, float phiStep)
{
  vtkSmartPointer<vtkPoints> points =
    vtkSmartPointer<vtkPoints>::New();

  vtkSmartPointer<vtkUnsignedCharArray> colors =
    vtkSmartPointer<vtkUnsignedCharArray>::New();
  colors->SetNumberOfComponents(3);
  colors->SetName("Colors");

  itk::ImageRegionConstIteratorWithIndex<RGBDIImageType> imageIterator(rgbdiImage,rgbdiImage->GetLargestPossibleRegion());

  while(!imageIterator.IsAtEnd())
    {
    // Get the value of the current pixel
    RGBDIImageType::PixelType pixel = imageIterator.Get();
    unsigned char rgb[3];
    rgb[0] = pixel[0];
    rgb[1] = pixel[1];
    rgb[2] = pixel[2];
    colors->InsertNextTupleValue(rgb);

    float theta = thetaStart + imageIterator.GetIndex()[0] * thetaStep;
    float phi = phiStart + imageIterator.GetIndex()[1] * phiStep;
    float rho = pixel[3];

    float x = rho * sin(theta) * cos(phi);
    float y = rho * sin(theta) * sin(phi);
    float z = rho * sin(phi);

    points->InsertNextPoint(x,y,z);

    ++imageIterator;
    }

  vtkSmartPointer<vtkPolyData> polydata =
    vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints(points);
  polydata->GetPointData()->SetScalars(colors);

  vtkSmartPointer<vtkVertexGlyphFilter> vertexGlyphFilter =
    vtkSmartPointer<vtkVertexGlyphFilter>::New();
  vertexGlyphFilter->AddInput(polydata);
  vertexGlyphFilter->Update();

  pointCloud->ShallowCopy(vertexGlyphFilter->GetOutput());
}