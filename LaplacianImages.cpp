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
  
  ptxImage.WriteXLaplacian(outputPrefix);
  ptxImage.WriteYLaplacian(outputPrefix);
  ptxImage.WriteZLaplacian(outputPrefix);
  ptxImage.WriteXYZLaplacian(outputPrefix);
  return EXIT_SUCCESS;
}
