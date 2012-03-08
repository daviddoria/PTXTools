#include "PTXImage.h"
#include "PTXReader.h"

int main (int argc, char *argv[])
{
  if(argc != 3)
    {
    std::cout << "Required arguments: InputFilename(ptx) OutputFilename(png)" << std::endl;
    return EXIT_FAILURE;
    }

  std::string inputFilename = argv[1];
  std::string outputFilename = argv[2];

  PTXImage ptxImage = PTXReader::Read(inputFilename);
  ptxImage.WriteIntensityImage(outputFilename);

  return EXIT_SUCCESS;
}
