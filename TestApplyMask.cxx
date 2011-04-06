#include "PTXImage.h"

#include "itkImage.h"
#include "itkImageFileReader.h"

int main (int argc, char *argv[])
{
  if(argc != 4)
    {
    std::cout << "Required arguments: InputFilename(ptx) MaskFilename(png) OutputFilename(png)" << std::endl;
    return EXIT_FAILURE;
    }
  
  std::string inputFilename = argv[1];
  std::string maskFilename = argv[2];
  std::string outputFilename = argv[3];

  typedef itk::ImageFileReader<itk::Image<unsigned char, 2> > ReaderType;
  ReaderType::Pointer maskReader = ReaderType::New();
  maskReader->SetFileName(maskFilename);
  maskReader->Update();
  
  PTXImage ptxImage;
  ptxImage.ReadFile(inputFilename);
  ptxImage.ApplyMask(maskReader->GetOutput());

  ptxImage.WriteRGBImage(outputFilename);
  ptxImage.WritePointCloud("pointcloud.vtp");

  return EXIT_SUCCESS;
}
