#include "PTXImage.h"
#include "PTXReader.h"

#include "itkImageFileReader.h"

int main (int argc, char *argv[])
{
  if(argc != 4)
    {
    std::cout << "Required arguments: InputFilename(ptx) NewDepthImage(mhd) OutputFilename(vtp)" << std::endl;
    return EXIT_FAILURE;
    }

  std::string inputFilename = argv[1];
  std::string newDepthImageFilename = argv[2];
  std::string outputFilename = argv[3];

  PTXImage ptxImage = PTXReader::Read(inputFilename);

  typedef itk::Image<itk::CovariantVector<float, 3>, 2> ImageType;
  typedef itk::ImageFileReader<ImageType> ImageReaderType;
  ImageReaderType::Pointer reader = ImageReaderType::New();
  reader->SetFileName(newDepthImageFilename);
  reader->Update();

  ptxImage.ReplaceXYZ(reader->GetOutput());

  ptxImage.WritePTX(outputFilename);

  return EXIT_SUCCESS;
}
