#ifndef PANE3D_H
#define PANE3D_H

#include "Pane.h"

#include <vtkPolyData.h>

struct Pane3D : public Pane
{
  Pane3D(QVTKWidget* qvtkWidget);
  
  vtkSmartPointer<vtkActor> PointCloudActor;
  vtkSmartPointer<vtkPolyDataMapper> PointCloudMapper;
  vtkSmartPointer<vtkPolyData> PointCloud;
};

#endif
