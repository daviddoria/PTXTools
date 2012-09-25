#ifndef PTXPIXEL_H
#define PTXPIXEL_H

#include "itkNumericTraits.h"

#include <iostream>

/** This class stores properties of a point in a PTX grid, including
    the color of the point (R,G,B), the Intensity of the return, and the coordinate (X,Y,Z) of the point.
    It also contains a flag, Valid, to track if the point is a valid return or not.*/
struct PTXPixel
{
  /** Constructor. */
  PTXPixel();

  /** A flag to track if the point is valid. */
  bool Valid;

  /** The color of the point. */
  unsigned char R;
  unsigned char G;
  unsigned char B;

  /** The Intensity of the reflection at the point. */
  float Intensity;

  /** The Coordinate of the point. */
  float X;
  float Y;
  float Z;

  /** Get the Coordinate of the point. */
  float GetCoordinate(const unsigned int coordinate) const;

  /** Output operator for PTXPixel */
  friend std::ostream& operator<<(std::ostream& output,  const PTXPixel &pixel);

  /** Compute the "left/right" angle of the point */
  float GetTheta() const;

  /** Compute the "up/down" angle of the point */
  float GetPhi() const;

  /** Compute depth (assuming scanner is at (0,0,0) ) */
  float GetDepth() const;

  /** Determine if two points are identical. */
  bool operator==(const PTXPixel &pixel) const;

  /** Determine if two points are not identical. */
  bool operator!=(const PTXPixel &pixel) const;

  /** Check if all coordinates of a point are zero. */
  bool IsZero() const;

  /** Check if a point is valid and not zero. */
  bool IsValid() const;
};

// This is needed to allow the TileImagesFilter to work with PTXPixel images.
namespace itk
{
template<>
class NumericTraits< PTXPixel >
{
private:
public:

  typedef PTXPixel Self;

  static const Self ZeroValue()
  {
    PTXPixel temp;
    return temp;
  }

  static unsigned int GetLength(Self)
  {
    return 1;
  }

  typedef PTXPixel PrintType;
};
} // end namespace itk

#endif
