#ifndef PTXPIXEL_H
#define PTXPIXEL_H

#include <iostream>

struct PTXPixel
{
  // This struct contains the color of the point (R,G,B), the Intensity of the return, and the coordinate (X,Y,Z) of the point.
  // It also contains a flag, Valid, to track if the point is a valid return or not.

  PTXPixel();

  // Is the point valid?
  bool Valid;

  // Color
  unsigned char R;
  unsigned char G;
  unsigned char B;

  // Intensity
  float Intensity;

  // Coordinate
  float X;
  float Y;
  float Z;

  // Output operator for PTXPixel
  friend std::ostream& operator<<(std::ostream& output,  const PTXPixel &pixel);

  // Compute the "left/right" angle of the point
  float GetTheta();

  // Compute the "up/down" angle of the point
  float GetPhi();

  // Compute depth (assuming scanner is at (0,0,0) )
  float GetDepth();

};

#endif