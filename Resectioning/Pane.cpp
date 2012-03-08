#include "Pane.h"

#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkSmartPointer.h>

Pane::Pane(QVTKWidget* inputQVTKWidget)
{
  this->Renderer = vtkSmartPointer<vtkRenderer>::New();
  this->SelectionStyle = NULL;

  this->qvtkWidget = inputQVTKWidget;
  this->qvtkWidget->GetRenderWindow()->AddRenderer(this->Renderer);
}

Pane::~Pane()
{

}

void Pane::Refresh()
{
  this->qvtkWidget->GetRenderWindow()->Render();
}