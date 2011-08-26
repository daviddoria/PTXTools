#include "PTXImage.h"

int main (int argc, char *argv[])
{
  if(argc != 2)
    {
    std::cerr << "Required arguments: InputFilename.ptx" << std::endl;
    return EXIT_FAILURE;
    }

  std::string inputFileName = argv[1];

  std::cout << "Input filename: " << inputFileName << std::endl;

  PTXImage ptxImage;
  ptxImage.ReadFile(inputFileName);

  std::cout << "Width: " << ptxImage.GetWidth() << std::endl;
  std::cout << "Height: " << ptxImage.GetHeight() << std::endl;

  std::cout << "Valid points: " << ptxImage.CountValidPoints() << std::endl;
  
  return EXIT_SUCCESS;
}
