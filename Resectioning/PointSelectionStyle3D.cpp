/*=========================================================================
 *
 *  Copyright David Doria 2011 daviddoria@gmail.com
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

#include "PointSelectionStyle3D.h"

// VTK
#include <vtkAbstractPicker.h>
#include <vtkActor.h>
#include <vtkActor2D.h>
#include <vtkCamera.h>
#include <vtkFollower.h>
#include <vtkGlyph3D.h>
#include <vtkLabeledDataMapper.h>
#include <vtkObjectFactory.h>
#include <vtkPointPicker.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRendererCollection.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkTextProperty.h>
#include <vtkVectorText.h>

// STL
#include <sstream>
#include <stdexcept>

// Custom
#include "ResectioningHelpers.h"

// Submodules
#include "Helpers/Helpers.h"

vtkStandardNewMacro(PointSelectionStyle3D);

PointSelectionStyle3D::PointSelectionStyle3D()
{
  this->MarkerRadius = .05;

  // Create a sphere to use as the dot
  this->DotSource = vtkSmartPointer<vtkSphereSource>::New();
  this->DotSource->SetRadius(this->MarkerRadius);
  this->DotSource->Update();

  SelectedPoints = vtkSmartPointer<vtkPoints>::New();
  SelectedPointsPolyData = vtkSmartPointer<vtkPolyData>::New();
  SelectedPointsPolyData->SetPoints(SelectedPoints);

  Glyph3D = vtkSmartPointer<vtkGlyph3D>::New();
  Glyph3D->SetSourceData(DotSource->GetOutput());
  Glyph3D->SetInputData(SelectedPointsPolyData);

  SelectedPointsMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  SelectedPointsMapper->SetInputConnection(Glyph3D->GetOutputPort());
  SelectedPointsActor = vtkSmartPointer<vtkActor>::New();
  SelectedPointsActor->SetMapper(SelectedPointsMapper);


  // Create the text numbers
  LabeledDataMapper = vtkSmartPointer<vtkLabeledDataMapper>::New();
  LabeledDataMapper->SetInputData(SelectedPointsPolyData);
  LabeledDataMapper->GetLabelTextProperty()->SetFontSize(20);
  LabelActor = vtkSmartPointer<vtkActor2D>::New();
  LabelActor->SetMapper(LabeledDataMapper);

}

void PointSelectionStyle3D::Initialize()
{
  // This function must be called after calling selectionStyle->SetCurrentRenderer(renderer);
  if(!this->CurrentRenderer)
  {
    throw std::runtime_error("Initialize: CurrentRenderer is NULL!");
  }
  
  this->CurrentRenderer->AddViewProp(SelectedPointsActor);
  this->CurrentRenderer->AddViewProp(LabelActor);
}

void PointSelectionStyle3D::OnLeftButtonDown() 
{
  //std::cout << "Picking pixel: " << this->Interactor->GetEventPosition()[0] << " "
  //          << this->Interactor->GetEventPosition()[1] << std::endl;

  vtkPointPicker* picker = vtkPointPicker::SafeDownCast(this->Interactor->GetPicker());
  if(!picker)
  {
    throw std::runtime_error("PointSelectionStyle3D::OnLeftButtonDown() picker is NULL!");
  }
  
  int success = picker->Pick(this->Interactor->GetEventPosition()[0],
  //vtkPointPicker::SafeDownCast(this->Interactor->GetPicker())->Pick(this->Interactor->GetEventPosition()[0],
          this->Interactor->GetEventPosition()[1],
          0,  // always zero.
          this->CurrentRenderer);
  if(!success)
    {
    // Forward events
    vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
    return;
    }

  vtkPointPicker* pointPicker = vtkPointPicker::SafeDownCast(this->Interactor->GetPicker());
  if(!pointPicker)
    {
    throw std::runtime_error("pointPicker is invalid!");
    }

  if(pointPicker->GetDataSet() != this->Data)
    {
    std::cerr << "pointPicker->GetDataSet(): " << pointPicker->GetDataSet() << " this->Data: "
              << this->Data << std::endl;
    throw std::runtime_error("Did not pick from the correct data set!");
    }
  /*
  vtkIdType pointId = vtkPointPicker::SafeDownCast(this->Interactor->GetPicker())->GetPointId();
  double p[3];
  this->Data->GetPoint(pointId, p);
  
  //std::cout << "Picked point: " << pointId << std::endl;
  //std::cout << "Point: " << pointId << " should have coordinate: "
              << p[0] << " " << p[1] << " " << p[2] << std::endl;
  */
  double picked[3] = {0,0,0};
  
  pointPicker->GetPickPosition(picked);
  //std::cout << "Picked point with coordinate: " << picked[0] << " " << picked[1] << " " << picked[2] << std::endl;

  if(this->Interactor->GetShiftKey())
    {
    this->CurrentRenderer->GetActiveCamera()->SetFocalPoint(picked);
    }

  // Only select the point if control is held
  if(this->Interactor->GetControlKey())
    {
    AddNumber(picked);
    }

  // Forward events
  vtkInteractorStyleTrackballCamera::OnLeftButtonDown();

}

void PointSelectionStyle3D::AddNumber(const double p[3])
{
  // Create the text
  std::stringstream ss;
  ss << Coordinates.size();
  
  std::cout << "Added 3D keypoint: " << p[0] << " " << p[1] << " " << p[2] << std::endl;
  Coord3D coord;
  coord.x = p[0];
  coord.y = p[1];
  coord.z = p[2];
  Coordinates.push_back(coord);
  std::cout << "There are " << Coordinates.size() << " 3D coordinates." << std::endl;
  
  SelectedPoints->InsertNextPoint(p);
  SelectedPointsPolyData->SetPoints(SelectedPoints);
  SelectedPointsPolyData->Modified();
  
  std::cout << "There are " << SelectedPointsPolyData->GetNumberOfPoints() << " selected 3D points." << std::endl;
  for(vtkIdType i = 0; i < SelectedPointsPolyData->GetNumberOfPoints(); ++i)
  {
    double selectedPoint[3];
    SelectedPointsPolyData->GetPoint(i, selectedPoint);
    std::cout << "Selected point " << i << " : " << selectedPoint[0] << " "
              << selectedPoint[1] << " " << selectedPoint[2] << std::endl;
  }
  
  // Create the dots
  
  // Create the text
  LabeledDataMapper->Update();
  
}

void PointSelectionStyle3D::RemoveAll()
{
//   ResectioningHelpers::RemoveAllActors(this->Numbers, this->CurrentRenderer);
//   ResectioningHelpers::RemoveAllActors(this->Points, this->CurrentRenderer);

  this->Numbers.clear();
  this->Points.clear();
  this->Coordinates.clear();
}

void PointSelectionStyle3D::SetCurrentRenderer(vtkRenderer* const renderer)
{
  vtkInteractorStyle::SetCurrentRenderer(renderer);
}

void PointSelectionStyle3D::DeleteLastCorrespondence()
{
  this->CurrentRenderer->RemoveViewProp( this->Numbers[this->Numbers.size() - 1]);
  this->CurrentRenderer->RemoveViewProp( this->Points[this->Points.size() - 1]);
  this->Numbers.erase(this->Numbers.end()-1);
  this->Points.erase(this->Points.end()-1);
  this->Coordinates.erase(this->Coordinates.end()-1);
}

unsigned int PointSelectionStyle3D::GetNumberOfCorrespondences()
{
  return this->SelectedPoints->GetNumberOfPoints();
}

Coord3D PointSelectionStyle3D::GetCorrespondence(const unsigned int correspondenceId)
{
  double p[3];
  this->SelectedPoints->GetPoint(correspondenceId, p);
  Coord3D coord(p[0], p[1], p[2]);
  return coord;
}
