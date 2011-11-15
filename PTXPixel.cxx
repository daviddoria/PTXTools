#include "PTXPixel.h"

#include <cmath>

#include "itkPoint.h"

PTXPixel::PTXPixel() : Valid(false), X(0), Y(0), Z(0)
{

}

float PTXPixel::GetCoordinate(const unsigned int coordinate) const
{
  if(coordinate == 0)
    {
    return X;
    }
  else if(coordinate == 1)
    {
    return Y;
    }
  else if(coordinate == 2)
    {
    return Z;
    }
  else
    {
    std::cerr << "Invalid coordinate " << coordinate << " specified..." << std::endl;
    return 0;
    }
}

float PTXPixel::GetTheta() const
{
  // Compute the "left/right" angle
  return atan(X/Y);
}

float PTXPixel::GetPhi() const
{
  // Compute the "up/down" angle
  return atan(Z/sqrt(X*X + Y*Y));
}

float PTXPixel::GetDepth() const
{
  if(!this->Valid)
    {
    return 0.0;
    }

  itk::Point<float, 3> origin;
  origin.Fill(0);

  itk::Point<float, 3> point;
  point[0] = this->X;
  point[1] = this->Y;
  point[2] = this->Z;

  return origin.EuclideanDistanceTo(point);
}

std::ostream& operator<<(std::ostream& output, const PTXPixel &pixel)
{
  // Output all of the information about the point
  output << "Coordinate: " << pixel.X << " " << pixel.Y << " " << pixel.Z << std::endl;
  output << "Color: " << (int)pixel.R << " " << (int)pixel.G << " " << (int)pixel.B << std::endl;
  output << "Intensity: " << pixel.Intensity << std::endl << std::endl;
  return output;
}

bool PTXPixel::operator==(const PTXPixel &pixel) const
{
  if(this->R == pixel.R &&
    this->G == pixel.G &&
    this->B == pixel.B &&
    this->Valid == pixel.Valid &&
    this->Intensity == pixel.Intensity &&
    this->X == pixel.X &&
    this->Y == pixel.Y &&
    this->Z == pixel.Z)
    {
    return true;
    }
  return false;
}

bool PTXPixel::operator!=(const PTXPixel &pixel) const
{
  return !(*this == pixel);
}
