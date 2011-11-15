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

#include "PTXViewer.h"

#include "Helpers.h"

// ITK
#include <itkCastImageFilter.h>
#include <itkCovariantVector.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkLineIterator.h>
#include <itkNthElementImageAdaptor.h>

// VTK
#include <vtkCamera.h>
#include <vtkImageActor.h>
#include <vtkImageData.h>
#include <vtkInteractorStyleImage.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkSmartPointer.h>

// Qt
#include <QFileDialog>
#include <QLineEdit>
#include <QMessageBox>

#include <iostream>

MainWindow::MainWindow(QWidget *parent)
{
  // Setup the GUI and connect all of the signals and slots
  setupUi(this);

  this->BackgroundColor[0] = 0;
  this->BackgroundColor[1] = 0;
  this->BackgroundColor[2] = .5;

  this->CameraUp[0] = 0;
  this->CameraUp[1] = 1;
  this->CameraUp[2] = 0;

  // Add renderers - we flip the image by changing the camera view up because of the conflicting conventions used by ITK and VTK
  this->LeftRenderer = vtkSmartPointer<vtkRenderer>::New();
  this->LeftRenderer->GradientBackgroundOn();
  this->LeftRenderer->SetBackground(this->BackgroundColor);
  this->LeftRenderer->SetBackground2(1,1,1);
  this->LeftRenderer->GetActiveCamera()->SetViewUp(this->CameraUp);
  this->qvtkWidgetLeft->GetRenderWindow()->AddRenderer(this->LeftRenderer);

  this->RightRenderer = vtkSmartPointer<vtkRenderer>::New();
  this->RightRenderer->GradientBackgroundOn();
  this->RightRenderer->SetBackground(this->BackgroundColor);
  this->RightRenderer->SetBackground2(1,1,1);
  this->RightRenderer->GetActiveCamera()->SetViewUp(this->CameraUp);
  this->qvtkWidgetRight->GetRenderWindow()->AddRenderer(this->RightRenderer);

  // Setup interactor styles
  this->InteractorStyleImage = vtkSmartPointer<vtkInteractorStyleImage>::New();
  this->qvtkWidgetLeft->GetInteractor()->SetInteractorStyle(this->InteractorStyleImage);
  
  this->InteractorStyleTrackballCamera = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
  this->qvtkWidgetRight->GetInteractor()->SetInteractorStyle(this->InteractorStyleTrackballCamera);

  // Things for the 2D window
  this->ColorImage = PTXImage::RGBImageType::New();
  this->ColorImageData = vtkSmartPointer<vtkImageData>::New();
  this->ColorImageSlice = vtkSmartPointer<vtkImageSlice>::New();
  this->ColorImageSliceMapper = vtkSmartPointer<vtkImageSliceMapper>::New();
  this->ColorImageSliceMapper->SetInputConnection(this->ColorImageData->GetProducerPort());
  this->ColorImageSlice->SetMapper(this->ColorImageSliceMapper);
  this->ColorImageSlice->VisibilityOff();
  this->LeftRenderer->AddViewProp(this->ColorImageSlice);
    
  this->IntensityImage = PTXImage::FloatImageType::New();
  this->IntensityImageData = vtkSmartPointer<vtkImageData>::New();
  this->IntensityImageSlice = vtkSmartPointer<vtkImageSlice>::New();
  this->IntensityImageSliceMapper = vtkSmartPointer<vtkImageSliceMapper>::New();
  this->IntensityImageSliceMapper->SetInputConnection(this->IntensityImageData->GetProducerPort());
  this->IntensityImageSlice->SetMapper(this->IntensityImageSliceMapper);
  this->IntensityImageSlice->VisibilityOff();
  this->LeftRenderer->AddViewProp(this->IntensityImageSlice);
  
  this->DepthImage = PTXImage::FloatImageType::New();
  this->DepthImageData = vtkSmartPointer<vtkImageData>::New();
  this->DepthImageSlice = vtkSmartPointer<vtkImageSlice>::New();
  this->DepthImageSliceMapper = vtkSmartPointer<vtkImageSliceMapper>::New();
  this->DepthImageSliceMapper->SetInputConnection(this->DepthImageData->GetProducerPort());
  this->DepthImageSlice->SetMapper(this->DepthImageSliceMapper);
  this->DepthImageSlice->VisibilityOff();
  this->LeftRenderer->AddViewProp(this->DepthImageSlice);
  
  // Things for the 3D window
  this->PointsActor = vtkSmartPointer<vtkActor>::New();
  this->PointsPolyData = vtkSmartPointer<vtkPolyData>::New();
  this->PointsPolyDataMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  this->PointsPolyDataMapper->SetInputConnection(this->PointsPolyData->GetProducerPort());
  this->PointsActor->SetMapper(this->PointsPolyDataMapper);
  this->RightRenderer->AddViewProp(this->PointsActor);
  
  // Default GUI settings
  this->radRGB->setChecked(true);
}

void MainWindow::on_actionFlipImage_activated()
{
  this->CameraUp[1] *= -1;
  this->LeftRenderer->GetActiveCamera()->SetViewUp(this->CameraUp);
  this->RightRenderer->GetActiveCamera()->SetViewUp(this->CameraUp);
  this->Refresh();
}

void MainWindow::on_actionQuit_activated()
{
  exit(-1);
}

void MainWindow::on_actionOpenImage_activated()
{
  OpenFile();
}

void MainWindow::on_radRGB_clicked()
{
  Refresh();
}

void MainWindow::on_radDepth_clicked()
{
  Refresh();
}

void MainWindow::on_radIntensity_clicked()
{
  Refresh();
}

#if 0
void InnerWidget::actionFlip_Image_triggered()
{
  this->CameraUp[1] *= -1;
  this->LeftRenderer->GetActiveCamera()->SetViewUp(this->CameraUp);
  this->RightRenderer->GetActiveCamera()->SetViewUp(this->CameraUp);
  this->Refresh();
}
#endif

void MainWindow::OpenFile()
{
  std::cout << "OpenFile()" << std::endl;
  
  // Get a filename to open
  QString filename = QFileDialog::getOpenFileName(this,
     "Open PTX File", ".", "Image Files (*.ptx)");

  if(filename.isEmpty())
    {
    return;
    }

  // Read file
  this->PTX.ReadFile(filename.toStdString());

  // Convert the images into a VTK images for display.
  this->PTX.CreateRGBImage(this->ColorImage);
  Helpers::ITKRGBImageToVTKImage(this->ColorImage, this->ColorImageData);
  
  this->PTX.CreateDepthImage(this->DepthImage);
  Helpers::ITKScalarImageToScaledVTKImage<PTXImage::FloatImageType>(this->DepthImage, this->DepthImageData);
  
  this->PTX.CreateIntensityImage(this->IntensityImage);
  Helpers::ITKScalarImageToScaledVTKImage<PTXImage::FloatImageType>(this->IntensityImage, this->IntensityImageData);

  this->PTX.CreatePointCloud(this->PointsPolyData);
  
  Refresh();
  
  this->LeftRenderer->ResetCamera();
  this->RightRenderer->ResetCamera();
}

void MainWindow::Refresh()
{
  this->DepthImageSlice->SetVisibility(this->radDepth->isChecked());
  this->IntensityImageSlice->SetVisibility(this->radIntensity->isChecked());
  this->ColorImageSlice->SetVisibility(this->radRGB->isChecked());
  
  this->LeftRenderer->Render();
  this->RightRenderer->Render();
  this->qvtkWidgetRight->GetRenderWindow()->Render();
  this->qvtkWidgetLeft->GetRenderWindow()->Render();
  this->qvtkWidgetRight->GetInteractor()->Render();
  this->qvtkWidgetLeft->GetInteractor()->Render();
}
