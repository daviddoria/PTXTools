/*=========================================================================
 *
 *  Copyright David Doria 2012 daviddoria@gmail.com
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

#include "PTXReader.h"

// STL
#include <iostream>

// ITK
#include "itkIndex.h"

int main(int argc, char*argv[])
{
  if(argc != 2)
  {
    throw std::runtime_error("Required arguments: file.ptx");
  }

  std::string ptxFileName = argv[1];

  PTXImage ptxImage = PTXReader::Read(ptxFileName);

  // Bottom left of the scan
  itk::Index<2> origin = {{0,0}};

  // Top right of the scan
  itk::Index<2> farCorner = {{static_cast<itk::Index<2>::IndexValueType>(ptxImage.GetFullImage()->GetLargestPossibleRegion().GetSize()[0] - 1),
                              static_cast<itk::Index<2>::IndexValueType>(ptxImage.GetFullImage()->GetLargestPossibleRegion().GetSize()[1] - 1)}};
  // Origin
  float originPhi = ptxImage.GetPhi(origin);
  std::cout << "Origin phi: " << originPhi << std::endl;
  std::cout << "Origin theta: " << ptxImage.GetTheta(origin) << std::endl;
  std::cout << "Origin z coordinate: " << ptxImage.GetFullImage()->GetPixel(origin).Z << std::endl;
  std::cout << "Origin y coordinate: " << ptxImage.GetFullImage()->GetPixel(origin).Y << std::endl;

  // Far corner
  std::cout << "farCorner phi: " << ptxImage.GetPhi(farCorner) << std::endl;
  std::cout << "farCorner theta: " << ptxImage.GetTheta(farCorner) << std::endl;
  std::cout << "farCorner z coordinate: " << ptxImage.GetFullImage()->GetPixel(farCorner).Z << std::endl;
  std::cout << "farCorner y coordinate: " << ptxImage.GetFullImage()->GetPixel(farCorner).Y << std::endl;

  return 0;
}
