#include "PTXImage.h"

int main (int argc, char *argv[])
{
  if(argc != 3)
    {
    std::cerr << "Required arguments: InputFilename(ptx) OutputFilename(ptx)" << std::endl;
    return EXIT_FAILURE;
    }

  std::string inputFilename = argv[1];
  std::string outputFilename = argv[2];

  PTXImage ptxImage;
  ptxImage.ReadFile(inputFilename);
  PTXImage downsampled = ptxImage.Downsample(2);

  return EXIT_SUCCESS;
}
