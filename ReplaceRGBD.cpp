#include "PTXImage.h"

#include "itkImageFileReader.h"

int main (int argc, char *argv[])
{
  if(argc != 4)
    {
    std::cout << "Required arguments: InputFilename(ptx) NewRGBD(mhd) OutputFilename(vtp)" << std::endl;
    return EXIT_FAILURE;
    }

  std::string inputFilename = argv[1];
  std::string newRGBDFilename = argv[2];
  std::string outputFilename = argv[3];

  PTXImage ptxImage;
  ptxImage.ReadFile(inputFilename);

  typedef itk::Image<itk::CovariantVector<float, 4>, 2> RGBDImageType;
  typedef itk::ImageFileReader<RGBDImageType> RGBDImageReaderType;
  RGBDImageReaderType::Pointer reader = RGBDImageReaderType::New();
  reader->SetFileName(newRGBDFilename);
  reader->Update();

  ptxImage.ReplaceRGBD(reader->GetOutput());

  ptxImage.WritePointCloud(outputFilename);

  return EXIT_SUCCESS;
}
