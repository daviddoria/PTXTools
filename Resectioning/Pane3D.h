#ifndef PANE3D_H
#define PANE3D_H

#include "Pane.h"

#include <vtkPolyData.h>

class vtkInteractorStyleTrackballCamera;

struct Pane3D : public Pane
{
  Pane3D(QVTKWidget* qvtkWidget);
  
  vtkSmartPointer<vtkActor> PointCloudActor;
  vtkSmartPointer<vtkPolyDataMapper> PointCloudMapper;
  vtkSmartPointer<vtkPolyData> PointCloud;

  void SetPolyData(vtkPolyData* const polyData);

  vtkInteractorStyleTrackballCamera* InteractorStyle;
  
private:
  
};

#endif
