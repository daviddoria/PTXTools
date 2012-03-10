// STL
#include <fstream>
#include <iostream>
#include <vector>

// Eigen
#include <Eigen/Geometry>
#include <Eigen/Dense> // for Vector

// VTK
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>

// ITK
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageRegionConstIterator.h"

// Custom
#include "PTXImage.h"
#include "PTXReader.h"
#include "Resectioning.h"

int main(int argc, char *argv[])
{
  if(argc < 5)
    {
    std::cerr << "Required arguments: cameraMatrixFileName.txt ptxFileName.ptx imageFileName.png outputFileName.ptx" << std::endl;
    return EXIT_FAILURE;
    }
  // Parse arguments
  std::string cameraMatrixFileName = argv[1];
  std::string ptxFileName = argv[2];
  std::string imageFileName = argv[3];
  std::string outputFileName = argv[4];

  // Output arguments
  std::cout << "cameraMatrixFileName : " << cameraMatrixFileName << std::endl;
  std::cout << "ptxFileName : " << ptxFileName << std::endl;
  std::cout << "imageFileName : " << imageFileName << std::endl;
  std::cout << "outputFileName : " << outputFileName << std::endl;

  // Read the PTX file
  PTXImage ptxImage = PTXReader::Read(ptxFileName);

  // Read the camera matrix relating the input image to the PTX/scan/LiDAR file
  Eigen::MatrixXd P = Resectioning::ReadP(cameraMatrixFileName);
  std::cout << "P: " << P << std::endl;

  // Read the image
  typedef itk::ImageFileReader<PTXImage::RGBImageType> ReaderType;
  ReaderType::Pointer imageReader = ReaderType::New();
  imageReader->SetFileName(imageFileName);
  imageReader->Update();

  Resectioning::ResectionNaive(P, ptxImage, imageReader->GetOutput());

  FilePrefix ptxPrefix(outputFileName);
  ptxImage.WritePTX(ptxPrefix);

  return EXIT_SUCCESS;
}
