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

/* This is the main GUI class of this project. It is a QMainWindow
 * so that we can use a File menu. It contains an instance of our main functional
 * class ImageGraphCutBase and our custom scribble interactor style vtkGraphCutInteractorStyle.
 * It also contains a CProgressThread so that we can display a progress bar in marquee
 * mode during long computations.
*/

#ifndef PTXViewerWidget_H
#define PTXViewerWidget_H

// Qt
#include "ui_PTXViewerWidget.h"

// VTK
#include <vtkActor.h>
#include <vtkImageActor.h>
#include <vtkImageData.h>
#include <vtkImageSlice.h>
#include <vtkImageSliceMapper.h>
#include <vtkInteractorStyleImage.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkSmartPointer.h>

// Custom
#include "ImageLayer.h"
#include "PTXImage.h"

class PTXViewerWidget : public QMainWindow, private Ui::PTXViewerWidget
{
Q_OBJECT
public:
  PTXViewerWidget(QWidget *parent = 0);

public slots:
  // Menu items
  // File menu
  void on_actionOpenImage_activated();
  void on_actionQuit_activated();
  void on_actionFlipImage_activated();
  
  // Export menu
  void on_actionExportRGBImage_activated();
  void on_actionExportRGBDImage_activated();
  void on_actionExportRGBDVImage_activated();
  void on_actionExportIntensityImage_activated();
  void on_actionExportDepthImage_activated();
  void on_actionExportValidityImage_activated();

  // Image display radio buttons.
  void on_radRGB_clicked();
  void on_radDepth_clicked();
  void on_radIntensity_clicked();
  void on_radValidity_clicked();

  // Use a QFileDialog to get a filename, then open the specified file as a greyscale or color image, depending on which type the user has specified through the file menu.
  void OpenFile();
  
protected:

  // Things for the 2D window
  vtkSmartPointer<vtkInteractorStyleImage> InteractorStyleImage;
  vtkSmartPointer<vtkRenderer> LeftRenderer;
  
  ImageLayer<PTXImage::RGBImageType> ColorImageLayer;
  
  ImageLayer<PTXImage::FloatImageType> IntensityImageLayer;
  
  ImageLayer<PTXImage::FloatImageType> DepthImageLayer;
  
  ImageLayer<PTXImage::UnsignedCharImageType> ValidityImageLayer;
  
  // Things for the 3D window
  vtkSmartPointer<vtkInteractorStyleTrackballCamera> InteractorStyleTrackballCamera;
  vtkSmartPointer<vtkRenderer> RightRenderer;
  vtkSmartPointer<vtkPolyData> PointsPolyData;
  vtkSmartPointer<vtkPolyDataMapper> PointsPolyDataMapper;
  vtkSmartPointer<vtkActor> PointsActor;
  
  // Refresh both renderers and render windows
  void Refresh();
  
  // Allows the background color to be changed
  double BackgroundColor[3];

  // Allows the image to be flipped so that it is "right side up"
  double CameraUp[3];

  // Store the PTX file when it is opened.
  PTXImage PTX;
  
  template <typename TImage>
  void SaveImage(const typename TImage::Pointer image);
};

#include "PTXViewerWidget.hxx"

#endif
