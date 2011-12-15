#include "PTXImage.h"

int main (int argc, char *argv[])
{
  // Verify arguments
  if(argc != 4)
    {
    std::cerr << "Required arguments: InputFilename(ptx) DownsampleFactor OutputFilename(ptx)" << std::endl;
    return EXIT_FAILURE;
    }

  // Parse arguments
  std::string inputFilename = argv[1];

  std::stringstream ss;
  ss << argv[2];
  unsigned int factor;
  ss >> factor;

  std::string outputFilename = argv[3];

  // Output arguments
  std::cout << "Input: " << inputFilename << std::endl;
  std::cout << "Factor: " << factor << std::endl;
  std::cout << "Output: " << outputFilename << std::endl;

  PTXImage ptxImage;
  ptxImage.ReadFile(inputFilename);
  PTXImage downsampled = ptxImage.Downsample(factor);

  downsampled.WritePTX(outputFilename);

  return EXIT_SUCCESS;
}
