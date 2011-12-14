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

#include "PTXViewerWidget.h"

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

PTXViewerWidget::PTXViewerWidget(QWidget *parent)
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
  this->LeftRenderer->AddViewProp(this->ColorImageLayer.ImageSlice);
  this->LeftRenderer->AddViewProp(this->DepthImageLayer.ImageSlice);
  this->LeftRenderer->AddViewProp(this->ValidityImageLayer.ImageSlice);
  this->LeftRenderer->AddViewProp(this->IntensityImageLayer.ImageSlice);

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

void PTXViewerWidget::on_actionFlipImage_activated()
{
  this->CameraUp[1] *= -1;
  this->LeftRenderer->GetActiveCamera()->SetViewUp(this->CameraUp);
  this->RightRenderer->GetActiveCamera()->SetViewUp(this->CameraUp);
  this->Refresh();
}

// File menu
void PTXViewerWidget::on_actionQuit_activated()
{
  exit(-1);
}

void PTXViewerWidget::on_actionOpenImage_activated()
{
  OpenFile();
}

// Image display radio buttons.
void PTXViewerWidget::on_radRGB_clicked()
{
  Refresh();
}

void PTXViewerWidget::on_radDepth_clicked()
{
  Refresh();
}

void PTXViewerWidget::on_radIntensity_clicked()
{
  Refresh();
}

void PTXViewerWidget::on_radValidity_clicked()
{
  Refresh();
}

// Export menu
void PTXViewerWidget::on_actionExportRGBImage_activated()
{
  SaveImage<PTXImage::RGBImageType>(this->ColorImageLayer.Image);
}

void PTXViewerWidget::on_actionExportRGBDImage_activated()
{
  PTXImage::RGBDImageType::Pointer rgbdImage = PTXImage::RGBDImageType::New();
  this->PTX.CreateRGBDImage(rgbdImage);
  SaveImage<PTXImage::RGBDImageType>(rgbdImage);
}

void PTXViewerWidget::on_actionExportRGBDVImage_activated()
{
  PTXImage::RGBDVImageType::Pointer rgbdvImage = PTXImage::RGBDVImageType::New();
  this->PTX.CreateRGBDVImage(rgbdvImage);
  SaveImage<PTXImage::RGBDVImageType>(rgbdvImage);
}

void PTXViewerWidget::on_actionExportIntensityImage_activated()
{
  SaveImage<PTXImage::FloatImageType>(this->IntensityImageLayer.Image);
}

void PTXViewerWidget::on_actionExportDepthImage_activated()
{
  SaveImage<PTXImage::FloatImageType>(this->DepthImageLayer.Image);
}

void PTXViewerWidget::on_actionExportValidityImage_activated()
{
  SaveImage<PTXImage::MaskImageType>(this->ValidityImageLayer.Image);
}

void PTXViewerWidget::on_actionExportPointCloud_activated()
{
  QString fileName = QFileDialog::getSaveFileName(this, "Save File", "", "VTP Files (*.vtp)");

  //DebugMessage<std::string>("Got filename: ", fileName.toStdString());
  if(fileName.toStdString().empty())
    {
    std::cout << "Filename was empty." << std::endl;
    return;
    }
    
  this->PTX.WritePointCloud(fileName.toStdString());
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

void PTXViewerWidget::OpenFile()
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
  this->PTX.CreateRGBImage(this->ColorImageLayer.Image);
  Helpers::ITKRGBImageToVTKImage(this->ColorImageLayer.Image, this->ColorImageLayer.ImageData);
  
  this->PTX.CreateDepthImage(this->DepthImageLayer.Image);
  Helpers::ITKScalarImageToScaledVTKImage<PTXImage::FloatImageType>(this->DepthImageLayer.Image, this->DepthImageLayer.ImageData);
  
  this->PTX.CreateIntensityImage(this->IntensityImageLayer.Image);
  Helpers::ITKScalarImageToScaledVTKImage<PTXImage::FloatImageType>(this->IntensityImageLayer.Image, this->IntensityImageLayer.ImageData);
  
  this->PTX.CreateValidityImage(this->ValidityImageLayer.Image);
  Helpers::ITKScalarImageToScaledVTKImage<PTXImage::UnsignedCharImageType>(this->ValidityImageLayer.Image, this->ValidityImageLayer.ImageData);

  this->PTX.CreatePointCloud(this->PointsPolyData);
  
  Refresh();
  
  this->LeftRenderer->ResetCamera();
  this->RightRenderer->ResetCamera();
}

void PTXViewerWidget::Refresh()
{
  this->DepthImageLayer.ImageSlice->SetVisibility(this->radDepth->isChecked());
  this->IntensityImageLayer.ImageSlice->SetVisibility(this->radIntensity->isChecked());
  this->ColorImageLayer.ImageSlice->SetVisibility(this->radRGB->isChecked());
  this->ValidityImageLayer.ImageSlice->SetVisibility(this->radValidity->isChecked());
  
  this->LeftRenderer->Render();
  this->RightRenderer->Render();
  this->qvtkWidgetRight->GetRenderWindow()->Render();
  this->qvtkWidgetLeft->GetRenderWindow()->Render();
  this->qvtkWidgetRight->GetInteractor()->Render();
  this->qvtkWidgetLeft->GetInteractor()->Render();
}
