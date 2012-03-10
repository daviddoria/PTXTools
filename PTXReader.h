/*
Copyright (C) 2011 David Doria, daviddoria@gmail.com

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

#ifndef PTXREADER_H
#define PTXREADER_H

// Custom
#include "PTXImage.h"

class PTXReader
{
public:
  static PTXImage Read(const std::string& filename);

  void SetFileName(const std::string& filename);

  void Read();

  PTXImage GetOutput();

private:
  std::string FileName;
  PTXImage ptxImage;
  
};

#endif
