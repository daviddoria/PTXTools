#include "Pane2D.h"

#include <vtkImageData.h>
#include <vtkImageProperty.h>
#include <vtkImageSliceMapper.h>
#include <vtkImageSlice.h>
#include <vtkRenderWindow.h>

#include "PointSelectionStyle2D.h"

Pane2D::Pane2D(QVTKWidget* qvtkWidget) : Pane(qvtkWidget)
{
  this->ImageData = vtkSmartPointer<vtkImageData>::New();

  this->ImageSliceMapper = vtkSmartPointer<vtkImageSliceMapper>::New();
  
  this->ImageSlice = vtkSmartPointer<vtkImageSlice>::New();
  // Make the pixels sharp instead of blurry when zoomed
  this->ImageSlice->GetProperty()->SetInterpolationTypeToNearest();
  
  this->CameraLeftToRightVector.resize(3);
  this->CameraLeftToRightVector[0] = -1;
  this->CameraLeftToRightVector[1] = 0;
  this->CameraLeftToRightVector[2] = 0;

  this->CameraBottomToTopVector.resize(3);
  this->CameraBottomToTopVector[0] = 0;
  this->CameraBottomToTopVector[1] = 1;
  this->CameraBottomToTopVector[2] = 0;

  this->InteractorStyle = NULL;
}

void Pane2D::SetCameraPosition()
{
  double leftToRight[3] = {this->CameraLeftToRightVector[0], this->CameraLeftToRightVector[1], this->CameraLeftToRightVector[2]};
  double bottomToTop[3] = {this->CameraBottomToTopVector[0], this->CameraBottomToTopVector[1], this->CameraBottomToTopVector[2]};

  static_cast<PointSelectionStyle2D*>(this->InteractorStyle)->SetImageOrientation(leftToRight, bottomToTop);
  static_cast<PointSelectionStyle2D*>(this->InteractorStyle)->GetCurrentRenderer()->ResetCamera();
  static_cast<PointSelectionStyle2D*>(this->InteractorStyle)->GetCurrentRenderer()->ResetCameraClippingRange();

  this->Renderer->ResetCamera();
  this->Renderer->ResetCameraClippingRange();
  this->Renderer->GetRenderWindow()->Render();
}

void Pane2D::FlipVertically()
{
  this->CameraBottomToTopVector[1] *= -1;
  SetCameraPosition();
}

void Pane2D::FlipHorizontally()
{
  this->CameraLeftToRightVector[0] *= -1;
  SetCameraPosition();
}
