#include "Pane3D.h"

#include <vtkActor.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkSmartPointer.h>

Pane3D::Pane3D(QVTKWidget* qvtkWidget) : Pane(qvtkWidget)
{
  this->PointCloud = vtkSmartPointer<vtkPolyData>::New();
  
  this->PointCloudMapper = vtkSmartPointer<vtkPolyDataMapper>::New();

  this->PointCloudActor = vtkSmartPointer<vtkActor>::New();
  this->PointCloudActor->SetMapper(this->PointCloudMapper);
  
  this->PointCloudMapper->SetInputData(this->PointCloud);

  this->Renderer->AddActor(this->PointCloudActor);

  this->InteractorStyle = NULL;
}

void Pane3D::SetPolyData(vtkPolyData* const polyData)
{
  this->PointCloud->DeepCopy(polyData);
}
