#include "PTXImage.h"
#include "PTXReader.h"

#include "itkImage.h"
#include "itkImageFileReader.h"

int main (int argc, char *argv[])
{
  if(argc != 4)
    {
    std::cout << "Required arguments: InputFilename(ptx) MaskFilename(png) OutputFilePrefix" << std::endl;
    return EXIT_FAILURE;
    }

  std::string inputFilename = argv[1];
  std::string maskFilename = argv[2];
  std::string outputFilePrefix = argv[3];

  typedef itk::ImageFileReader<itk::Image<unsigned char, 2> > ReaderType;
  ReaderType::Pointer maskReader = ReaderType::New();
  maskReader->SetFileName(maskFilename);
  maskReader->Update();

  PTXImage ptxImage = PTXReader::Read(inputFilename);
  ptxImage.ApplyMask(maskReader->GetOutput());

  ptxImage.WriteRGBImage(outputFilePrefix);

  std::stringstream pointsss;
  pointsss << outputFilePrefix << ".vtp";
  ptxImage.WritePointCloud(pointsss.str());

  return EXIT_SUCCESS;
}
