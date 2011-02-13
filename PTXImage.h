#ifndef PTXImage_h
#define PTXImage_h
#include <boost/graph/graph_concepts.hpp>

#include "itkImage.h"

struct PTXPixel
{
  unsigned char R;
  unsigned char G;
  unsigned char B;
  float Intensity;
  float X;
  float Y;
  float Z;
  friend std::ostream& operator<<(std::ostream& output,  const PTXPixel &pixel);
  float GetTheta()
  {
    return atan(X/Y);
  }

  float GetPhi()
  {
    return atan(Z/sqrt(X*X + Y*Y));
  }
};

std::ostream& operator<<(std::ostream& output, const PTXPixel &pixel)
{
  output << "Coordinate: " << pixel.X << " " << pixel.Y << " " << pixel.Z << std::endl;
  output << "Color: " << (int)pixel.R << " " << (int)pixel.G << " " << (int)pixel.B << std::endl;
  output << "Intensity: " << pixel.Intensity << std::endl << std::endl;
  return output;
}

typedef itk::Image<PTXPixel, 2> FullImageType;

class PTXImage
{
public:
  PTXImage();

  void WriteRGBImage(std::string filename);

  void WritePointCloud(std::string filename);

  void ReadFile(std::string filename);

  FullImageType::Pointer FullImage;

};

#include "PTXImage.txx"

#endif