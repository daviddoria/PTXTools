#ifndef PANE2D_H
#define PANE2D_H

#include "Pane.h"

#include "Types.h"

class vtkImageSliceMapper;
class vtkImageSlice;
class vtkInteractorStyleImage;

struct Pane2D : public Pane
{
  Pane2D(QVTKWidget* qvtkWidget);

  FloatVectorImageType::Pointer Image;
  vtkSmartPointer<vtkImageData> ImageData;

  vtkSmartPointer<vtkImageSliceMapper> ImageSliceMapper;
  vtkSmartPointer<vtkImageSlice> ImageSlice;

  void FlipVertically();
  void FlipHorizontally();

  vtkInteractorStyleImage* InteractorStyle;
  
private:
  std::vector<float> CameraLeftToRightVector;
  std::vector<float> CameraBottomToTopVector;
  void SetCameraPosition();

  
};

#endif
