#include "PTXImage.h"

#include "itkImageFileReader.h"

int main (int argc, char *argv[])
{
  if(argc != 4)
    {
    std::cout << "Required arguments: InputFilename(ptx) NewRGB(mhd) OutputFilename(ptx)" << std::endl;
    return EXIT_FAILURE;
    }

  std::string inputFilename = argv[1];
  std::string newDepthImageFilename = argv[2];
  std::string outputFilename = argv[3];

  PTXImage ptxImage;
  ptxImage.ReadFile(inputFilename);

  typedef itk::Image<itk::CovariantVector<float, 3>, 2> ImageType;
  typedef itk::ImageFileReader<ImageType> ImageReaderType;
  ImageReaderType::Pointer reader = ImageReaderType::New();
  reader->SetFileName(newDepthImageFilename);
  reader->Update();
  
  if(reader->GetOutput()->GetLargestPossibleRegion().GetSize() != ptxImage.GetSize())
    {
    std::cerr << "RGB image must be the same size as PTX image!" << std::endl;
    exit(-1);
    }

  ptxImage.ReplaceRGB(reader->GetOutput());

  ptxImage.WritePTX(outputFilename);

  return EXIT_SUCCESS;
}
