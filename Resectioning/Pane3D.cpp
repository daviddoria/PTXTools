#include "Pane3D.h"

#include <vtkActor.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkSmartPointer.h>

Pane3D::Pane3D(QVTKWidget* qvtkWidget) : Pane(qvtkWidget)
{
  this->PointCloudActor = vtkSmartPointer<vtkActor>::New();
  this->PointCloudMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  this->PointCloud = vtkSmartPointer<vtkPolyData>::New();
}
