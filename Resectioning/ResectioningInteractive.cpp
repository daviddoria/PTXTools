/*=========================================================================
 *
 *  Copyright David Doria 2011 daviddoria@gmail.com
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

#include <QApplication>
#include <QCleanlooksStyle>

#include <stdexcept>

#include "ResectioningWidget.h"

// Run with
// camera.png scan.ptx
// or
// camera.png 2D.txt scan.ptx 3D.txt

int main( int argc, char** argv )
{
  QApplication app( argc, argv );

  QApplication::setStyle(new QCleanlooksStyle);

  ResectioningWidget* resectioningWidget = NULL;

  if(argc == 1)
  {
    std::cout << "Using no arguments." << std::endl;
    resectioningWidget = new ResectioningWidget;
  }
  else if(argc == 3)
  {
    std::cout << "Using files arguments." << std::endl;
    std::string imageFileName = argv[1];
    std::string pointCloudFileName = argv[2];
    std::cout << "Image: " << imageFileName << std::endl;
    std::cout << "Point cloud: " << pointCloudFileName << std::endl;
    resectioningWidget = new ResectioningWidget(imageFileName, pointCloudFileName);
  }
  else if(argc == 5)
  {
    std::cout << "Using files+correspondences arguments." << std::endl;
    std::string imageFileName = argv[1];
    std::string imageCorrespondencesFile = argv[2];
    std::string pointCloudFileName = argv[3];
    std::string pointCloudCorrespondencesFile = argv[4];
    std::cout << "Image: " << imageFileName << std::endl;
    std::cout << "Image correspondences: " << imageCorrespondencesFile << std::endl;
    std::cout << "Point cloud: " << pointCloudFileName << std::endl;
    std::cout << "Point cloud correspondences: " << pointCloudCorrespondencesFile << std::endl;
    resectioningWidget = new ResectioningWidget(imageFileName, imageCorrespondencesFile,
                                                pointCloudFileName, pointCloudCorrespondencesFile);
  }
  else
  {
    throw std::runtime_error("Bad arguments!");
  }
  resectioningWidget->show();

  return app.exec();
}
