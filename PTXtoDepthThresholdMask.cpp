#include "PTXImage.h"
#include "PTXReader.h"

int main (int argc, char *argv[])
{
  if(argc != 4)
    {
    std::cout << "Required arguments: InputFilename(ptx) OutputFilename(png) depthThreshold" << std::endl;
    return EXIT_FAILURE;
    }

  std::string inputFilename = argv[1];
  std::string outputFilename = argv[2];
  std::string strDepthThreshold = argv[3];
  std::stringstream ss;
  ss << strDepthThreshold;
  float depthThreshold;
  ss >> depthThreshold;

  PTXImage ptxImage = PTXReader::Read(inputFilename);

  ptxImage.WriteDepthThresholdMask(outputFilename, depthThreshold);

  return EXIT_SUCCESS;
}
