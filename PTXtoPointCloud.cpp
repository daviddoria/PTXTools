#include "PTXImage.h"
#include "PTXReader.h"
#include "FilePrefix.h"

int main (int argc, char *argv[])
{
  if(argc != 3)
    {
    std::cout << "Required arguments: Input (PTX) OutputPrefix" << std::endl;
    return EXIT_FAILURE;
    }

  std::string inputFilename = argv[1];
  std::string outputPrefix = argv[2];

  FilePrefix prefix(outputPrefix);

  PTXImage ptxImage = PTXReader::Read(inputFilename);
  ptxImage.WritePointCloud(prefix);

  return EXIT_SUCCESS;
}
