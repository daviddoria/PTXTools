#ifndef Coord_H
#define Coord_H

struct Coord
{

};

struct Coord2D : public Coord
{
  float x,y;
  Coord2D() : x(0), y(0){}
  Coord2D(const float x_in, const float y_in) : x(x_in), y(y_in){}
};

struct Coord3D : public Coord
{
  float x,y,z;
  Coord3D() : x(0), y(0), z(0){}
  Coord3D(const float x_in, const float y_in, const float z_in) : x(x_in), y(y_in), z(z_in){}
};

#endif
