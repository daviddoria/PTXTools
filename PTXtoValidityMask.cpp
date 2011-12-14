#include "PTXImage.h"

int main (int argc, char *argv[])
{
  if(argc != 3)
    {
    std::cout << "Required arguments: InputFilename(ptx) OutputPrefix (output will be OutputPrefix_ValidityMask.png)" << std::endl;
    return EXIT_FAILURE;
    }

  std::string inputFilename = argv[1];
  std::string outputPrefix = argv[2];

  PTXImage ptxImage;
  ptxImage.ReadFile(inputFilename);
  ptxImage.WriteInvalidMask(outputPrefix);

  return EXIT_SUCCESS;
}
