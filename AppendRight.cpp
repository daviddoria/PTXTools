#include "PTXImage.h"

int main (int argc, char *argv[])
{
  if(argc != 4)
    {
    std::cerr << "Required arguments: InputFilename.ptx PTXtoAppend.ptx OutputPrefix" << std::endl;
    return EXIT_FAILURE;
    }

  std::string inputFilename = argv[1];
  std::string inputToAppendFilename = argv[2];
  std::string outputFileName = argv[3];

  std::cout << "Reading first image..." << std::endl;
  PTXImage ptxImage;
  ptxImage.ReadFile(inputFilename);

  std::cout << "Reading second image..." << std::endl;
  PTXImage ptxImageToAppend;
  ptxImageToAppend.ReadFile(inputToAppendFilename);

  std::cout << "Appending images..." << std::endl;
  ptxImage.AppendPTXRight(ptxImageToAppend);

  std::cout << "Writing output..." << std::endl;
  ptxImage.WritePTX(outputFileName);
  
  return EXIT_SUCCESS;
}
