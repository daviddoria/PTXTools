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

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

// Qt
#include "ui_PTXViewer.h"

// VTK
#include <vtkActor.h>
#include <vtkImageActor.h>
#include <vtkInteractorStyleImage.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkSmartPointer.h>

// Custom
#include "PTXImage.h"

class MainWindow : public QMainWindow, private Ui::MainWindow
{
Q_OBJECT
public:
  MainWindow(QWidget *parent = 0);

public slots:
  // Menu items
  void actionOpenImage_triggered();
  void actionQuit_triggered();
  void actionFlipImage_triggered();

  // Buttons, radio buttons, and sliders
  void on_radRGB_clicked();
  void on_radDepth_clicked();
  void on_radIntensity_clicked();

  // Use a QFileDialog to get a filename, then open the specified file as a greyscale or color image, depending on which type the user has specified through the file menu.
  void OpenFile();
  
  void DisplayRGB();
  void DisplayDepth();
  void DisplayIntensity();
  
protected:

  // Our scribble interactor style
  vtkSmartPointer<vtkInteractorStyleImage> InteractorStyleImage;
  vtkSmartPointer<vtkInteractorStyleTrackballCamera> InteractorStyleTrackballCamera;

  // The input and output image actors
  vtkSmartPointer<vtkImageActor> ImageActor;
  vtkSmartPointer<vtkActor> Actor;

  // The renderers
  vtkSmartPointer<vtkRenderer> LeftRenderer;
  vtkSmartPointer<vtkRenderer> RightRenderer;

  // Refresh both renderers and render windows
  void Refresh();

  // Allows the background color to be changed
  double BackgroundColor[3];

  // Allows the image to be flipped so that it is "right side up"
  double CameraUp[3];

  // Store the image when it is opened.
  PTXImage PTX;
};

#endif
