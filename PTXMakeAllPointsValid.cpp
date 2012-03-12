#include "PTXImage.h"
#include "PTXReader.h"

#include "itkImageFileReader.h"

int main (int argc, char *argv[])
{
  if(argc != 3)
    {
    std::cout << "Required arguments: InputFilename(ptx) OutputFilename(ptx)" << std::endl;
    return EXIT_FAILURE;
    }

  std::string inputFilename = argv[1];
  std::string outputFilename = argv[2];

  std::cout << "input: " << inputFilename << std::endl;
  std::cout << "output: " << outputFilename << std::endl;
  
  PTXImage ptxImage = PTXReader::Read(inputFilename);

  std::cout << "Input has " << ptxImage.CountInvalidPoints() << " invalid points." << std::endl;
  
  ptxImage.SetAllPointsToValid();

  std::cout << "Output has " << ptxImage.CountInvalidPoints() << " invalid points." << std::endl;
  
  ptxImage.WritePTX(outputFilename);

  return EXIT_SUCCESS;
}
