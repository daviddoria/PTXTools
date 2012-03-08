#ifndef PANE_H
#define PANE_H

#include "PointSelectionStyle.h"

#include <vtkSmartPointer.h>

#include <QVTKWidget.h>

class vtkRenderer;

struct Pane
{
  Pane(QVTKWidget* qvtkWidget);
  virtual ~Pane();
  
  vtkSmartPointer<vtkRenderer> Renderer;

  PointSelectionStyle* SelectionStyle;

  QVTKWidget* qvtkWidget;

  void Refresh();

  
};

#endif
