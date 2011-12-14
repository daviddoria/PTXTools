#include "PTXImage.h"

int main (int argc, char *argv[])
{
  if(argc != 7)
    {
    std::cerr << "Required arguments: InputFilename.ptx startX startY sizeX sizeY OutputPrefix" << std::endl;
    return EXIT_FAILURE;
    }

  std::string inputFileName = argv[1];
  std::string strStartX = argv[2];
  std::string strStartY = argv[3];
  std::string strSizeX = argv[4];
  std::string strSizeY = argv[5];
  std::string outputFileName = argv[6];

  unsigned int startX, startY, sizeX, sizeY;

  std::stringstream ss;
  ss.str(strStartX);
  ss >> startX;
  ss.clear();

  ss.str(strStartY);
  ss >> startY;
  ss.clear();

  ss.str(strSizeX);
  ss >> sizeX;
  ss.clear();

  ss.str(strSizeY);
  ss >> sizeY;
  ss.clear();

  std::cout << "Input filename: " << inputFileName << std::endl;
  std::cout << "Output filename: " << outputFileName << std::endl;
  std::cout << "start X: " << startX << std::endl;
  std::cout << "start Y: " << startY << std::endl;
  std::cout << "size X: " << sizeX << std::endl;
  std::cout << "size Y: " << sizeY << std::endl;
  
  itk::Index<2> start;
  start[0] = startX;
  start[1] = startY;

  itk::Size<2> size;
  size[0] = sizeX;
  size[1] = sizeY;

  itk::ImageRegion<2> region(start,size);

  std::cout << "Crop region: " << region << std::endl;

  std::cout << "Reading PTX..." << std::endl;
  PTXImage ptxImage;
  ptxImage.ReadFile(inputFileName);

  std::cout << "Cropping..." << std::endl;
  ptxImage.Crop(region);

  std::cout << "Writing output..." << std::endl;
  ptxImage.WritePTX(outputFileName);

  return EXIT_SUCCESS;
}
