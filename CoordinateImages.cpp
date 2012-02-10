#include "PTXImage.h"
#include "PTXReader.h"

#include "itkImageFileWriter.h"

int main (int argc, char *argv[])
{
  if(argc != 3)
    {
    std::cerr << "Required arguments: InputFilename(ptx) OutputPrefix" << std::endl;
    return EXIT_FAILURE;
    }

  std::string inputFilename = argv[1];
  std::string outputPrefix = argv[2];

  PTXImage ptxImage = PTXReader::Read(inputFilename);

  PTXImage::FloatImageType::Pointer xImage = ptxImage.GetXImage();
  PTXImage::FloatImageType::Pointer yImage = ptxImage.GetYImage();
  PTXImage::FloatImageType::Pointer zImage = ptxImage.GetZImage();

  typedef itk::ImageFileWriter< PTXImage::FloatImageType > SingleWriterType;

  std::stringstream ssX;
  ssX << outputPrefix << "_x.mhd";
  SingleWriterType::Pointer xWriter = SingleWriterType::New();
  xWriter->SetFileName(ssX.str());
  xWriter->SetInput(xImage);
  xWriter->Update();

  std::stringstream ssY;
  ssY << outputPrefix << "_y.mhd";
  SingleWriterType::Pointer yWriter = SingleWriterType::New();
  yWriter->SetFileName(ssY.str());
  yWriter->SetInput(yImage);
  yWriter->Update();

  std::stringstream ssZ;
  ssZ << outputPrefix << "_z.mhd";
  SingleWriterType::Pointer zWriter = SingleWriterType::New();
  zWriter->SetFileName(ssZ.str());
  zWriter->SetInput(zImage);
  zWriter->Update();

  // Write XYZ image
  typedef itk::ImageFileWriter< PTXImage::XYZImageType > FullWriterType;

  std::stringstream ssXYZ;
  ssXYZ << outputPrefix << "_xyz.mhd";
  FullWriterType::Pointer xyzWriter = FullWriterType::New();
  xyzWriter->SetFileName(ssXYZ.str());
  xyzWriter->SetInput(ptxImage.GetXYZImage());
  xyzWriter->Update();

  return EXIT_SUCCESS;
}
