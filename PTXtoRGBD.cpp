#include "PTXImage.h"

int main (int argc, char *argv[])
{
  if(argc != 3)
    {
    std::cout << "Required arguments: InputFilename.ptx OutputPrefix (will automaticall append _RGBD.mha)" << std::endl;
    return EXIT_FAILURE;
    }

  std::string inputFilename = argv[1];
  std::string outputPrefix = argv[2];

  PTXImage ptxImage;
  ptxImage.ReadFile(inputFilename);
  ptxImage.WriteRGBDImage(outputPrefix);

  return EXIT_SUCCESS;
}
