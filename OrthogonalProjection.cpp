#include "PTXImage.h"
#include "PTXReader.h"

#include <sstream>

int main (int argc, char *argv[])
{
  if(argc != 3)
    {
    std::cout << "Required arguments: InputFilename(ptx) OutputFilename(png)" << std::endl;
    return EXIT_FAILURE;
    }

  std::string inputFilename = argv[1];
  std::string outputPrefix = argv[2];

  PTXImage ptxImage = PTXReader::Read(inputFilename);

  //PTXImage::VectorType projectionAxis = ptxImage.GetPrincipalAxis();
  PTXImage::VectorType projectionAxis;
  projectionAxis[0] = .0036;
  projectionAxis[1] = -.0026;
  projectionAxis[2] = .9999;

  PTXImage orthoPTX = ptxImage.OrthogonalProjection(projectionAxis);

  orthoPTX.WriteEverything(outputPrefix);

  return EXIT_SUCCESS;
}
