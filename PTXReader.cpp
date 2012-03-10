/*
Copyright (C) 2010 David Doria, daviddoria@gmail.com

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "PTXReader.h"

#include <stdexcept>

void PTXReader::Read()
{
  this->ptxImage = Read(this->FileName);

}

PTXImage PTXReader::GetOutput()
{
  return this->ptxImage;
}

void PTXReader::SetFileName(const std::string& filename)
{
  this->FileName = filename;
}

PTXImage PTXReader::Read(const std::string& filename)
{
  PTXImage ptxImage;
  
  // Attempt to open the file
  std::ifstream infile;
  infile.open(filename.c_str());

  // Verify that the file was opened correctly
  if(!infile)
    {
    std::cout << "Could not open file " << filename << "!" << std::endl;
    throw std::runtime_error("File not found!");
    }

  // Read the header
  std::string line;

  unsigned int numberOfThetaPoints; // aka numberOfColumns
  unsigned int numberOfPhiPoints; // aka numberOfRows

  getline(infile, line);
  std::stringstream(line) >> numberOfThetaPoints;

  getline(infile, line);
  std::stringstream(line) >> numberOfPhiPoints;

  std::cout << "PhiPoints: " << numberOfPhiPoints << std::endl;
  std::cout << "ThetaPoints: " << numberOfThetaPoints << std::endl;

  // Setup the image
  itk::Size<2> size;
  size[0] = numberOfThetaPoints;
  size[1] = numberOfPhiPoints;

  itk::Index<2> start;
  start.Fill(0);

  itk::ImageRegion<2> region(start, size);

  ptxImage.FullImage->SetRegions(region);
  ptxImage.FullImage->Allocate();

  // Skip 8 lines (identity matrices)
  for(int i = 0; i < 8; i++)
    {
    getline(infile, line);
    }

  for(unsigned int theta = 0; theta < numberOfThetaPoints; theta++)
    {
    for(unsigned int phi = 0; phi < numberOfPhiPoints; phi++)
      {
      itk::Index<2> pixelIndex;
      pixelIndex[0] = theta;
      pixelIndex[1] = phi;

      getline(infile, line);

      float coordinate[3];
      float intensity;
      //unsigned char color[3];
      int color[3]; // must read these from the string stream as an int

      std::stringstream parsedLine(line);
      parsedLine >> coordinate[0] >> coordinate[1] >> coordinate[2] >> intensity
                 >> color[0] >> color[1] >> color[2];
      //std::cout << color[0] << " " << color[1] << " " << color[2] << std::endl;

      // Create the pixel with all of the properties
      PTXPixel pixel;
      pixel.X = coordinate[0];
      pixel.Y = coordinate[1];
      pixel.Z = coordinate[2];
      pixel.Intensity = intensity;
      pixel.R = color[0];
      pixel.G = color[1];
      pixel.B = color[2];

      // Check for a particular value (0 0 0 0.5 0 0 0) which indicates that the scanner did not receive a valid return.
      if(pixel.X == 0 && pixel.Y == 0 && pixel.Z == 0 &&
        intensity == 0.50 && pixel.R == 0 && pixel.G == 0 && pixel.B == 0)
        {
        pixel.Valid = false;
        }
      else
        {
        pixel.Valid = true;
        }

      // Set the pixel in the image
      ptxImage.FullImage->SetPixel(pixelIndex, pixel);
      //std::cout << "Pixel " << pixelIndex << " Valid? " << pixel.Valid << std::endl;
      }//end phi for loop
    }// end theta for loop

  // Close the input file
  infile.close();

  ptxImage.Backup();

  return ptxImage;
}
