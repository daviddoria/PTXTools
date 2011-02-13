#include "PTXImage.h"

int main (int argc, char *argv[])
{
  if(argc != 3)
    {
    std::cout << "Required arguments: Input (PTX) Output(vtp)" << std::endl;
    return EXIT_FAILURE;
    }

  std::string inputFilename = argv[1];
  std::string outputFilename = argv[2];

  PTXImage ptxImage;
  ptxImage.ReadFile(inputFilename);
  ptxImage.WritePointCloud(outputFilename);


  return EXIT_SUCCESS;
}
