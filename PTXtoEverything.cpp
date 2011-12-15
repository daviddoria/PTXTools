#include "PTXImage.h"

#include "itkImageFileWriter.h"

#include <sstream>

int main (int argc, char *argv[])
{
  if(argc != 3)
    {
    std::cout << "Required arguments: InputFilename(ptx) OutputPrefix" << std::endl;
    return EXIT_FAILURE;
    }

  std::string inputFilename = argv[1];
  std::string outputPrefix = argv[2];

  PTXImage ptxImage;
  ptxImage.ReadFile(inputFilename);

  // RGB
  ptxImage.WriteRGBImage(outputPrefix);

  // Point cloud
  ptxImage.WritePointCloud(outputPrefix);

  // Depth image
  ptxImage.WriteDepthImage(outputPrefix);

  // Depth image
  ptxImage.WriteIntensityImage(outputPrefix);

  ptxImage.WriteRGBDImage(outputPrefix);

  ptxImage.WriteRGBDIImage(outputPrefix);

  ptxImage.WriteInvalidMask(outputPrefix);

  ptxImage.WriteDepthThresholdMask(outputPrefix, 2.5);

  ptxImage.WriteX(outputPrefix);
  ptxImage.WriteY(outputPrefix);
  ptxImage.WriteZ(outputPrefix);

  ptxImage.WriteXYZ(outputPrefix);

  ptxImage.WriteXLaplacian(outputPrefix);
  ptxImage.WriteYLaplacian(outputPrefix);
  ptxImage.WriteZLaplacian(outputPrefix);

  ptxImage.WriteXYZLaplacian(outputPrefix);


  return EXIT_SUCCESS;
}
