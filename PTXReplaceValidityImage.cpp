#include "PTXImage.h"
#include "PTXReader.h"

#include "itkImageFileReader.h"

int main (int argc, char *argv[])
{
  if(argc != 4)
    {
    std::cerr << "Required arguments: InputFilename(ptx) ValidityImage(png) OutputFilename(ptx)" << std::endl;
    return EXIT_FAILURE;
    }

  std::string inputFilename = argv[1];
  std::string validityFilename = argv[2];
  std::string outputFilename = argv[3];

  PTXImage ptxImage = PTXReader::Read(inputFilename);

  typedef itk::Image<unsigned char, 2> MaskType;
  typedef itk::ImageFileReader<MaskType> MaskReaderType;
  MaskReaderType::Pointer reader = MaskReaderType::New();
  reader->SetFileName(validityFilename);
  reader->Update();
  
  ptxImage.ReplaceValidity(reader->GetOutput());

  ptxImage.WritePTX(outputFilename);
  return EXIT_SUCCESS;
}
