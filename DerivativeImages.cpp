#include "PTXImage.h"

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

  PTXImage ptxImage;
  ptxImage.ReadFile(inputFilename);
  
  PTXImage::FloatImageType::Pointer xDerivative = ptxImage.GetXDerivative();
  PTXImage::FloatImageType::Pointer yDerivative = ptxImage.GetYDerivative();
  PTXImage::FloatImageType::Pointer zDerivative = ptxImage.GetZDerivative();

  typedef  itk::ImageFileWriter< PTXImage::FloatImageType > WriterType;
  
  std::stringstream ssX;
  ssX << outputPrefix << "_x.mhd";
  WriterType::Pointer xWriter = WriterType::New();
  xWriter->SetFileName(ssX.str());
  xWriter->SetInput(xDerivative);
  xWriter->Update();
 
  std::stringstream ssY;
  ssY << outputPrefix << "_y.mhd";
  WriterType::Pointer yWriter = WriterType::New();
  yWriter->SetFileName(ssY.str());
  yWriter->SetInput(yDerivative);
  yWriter->Update();
  
  std::stringstream ssZ;
  ssZ << outputPrefix << "_z.mhd";
  WriterType::Pointer zWriter = WriterType::New();
  zWriter->SetFileName(ssZ.str());
  zWriter->SetInput(zDerivative);
  zWriter->Update();
  return EXIT_SUCCESS;
}
