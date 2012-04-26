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

#include "PointSelectionStyle2D.h"

// VTK
#include <vtkAbstractPicker.h>
#include <vtkActor2D.h>
#include <vtkCaptionActor2D.h>
#include <vtkCoordinate.h>
#include <vtkFollower.h>
#include <vtkObjectFactory.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataMapper2D.h>
#include <vtkProp.h>
#include <vtkProperty.h>
#include <vtkProperty2D.h>
#include <vtkRendererCollection.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkTextProperty.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkVectorText.h>

// STL
#include <sstream>

// Custom
#include "ResectioningHelpers.h"

// Submodules
#include "Helpers/Helpers.h"

vtkStandardNewMacro(PointSelectionStyle2D);
 
void PointSelectionStyle2D::OnLeftButtonDown() 
{
  //std::cout << "Picking pixel: " << this->Interactor->GetEventPosition()[0] << " " << this->Interactor->GetEventPosition()[1] << std::endl;
  this->Interactor->GetPicker()->Pick(this->Interactor->GetEventPosition()[0], 
		      this->Interactor->GetEventPosition()[1], 
		      0,  // always zero.
                      this->CurrentRenderer);

  double picked[3];
  this->Interactor->GetPicker()->GetPickPosition(picked);
  //std::cout << "Picked point with coordinate: " << picked[0] << " " << picked[1] << " " << picked[2] << std::endl;

  AddNumber(picked);
 
  // Forward events
  vtkInteractorStyleImage::OnLeftButtonDown();
}
/*
void PointSelectionStyle2D::RemoveAllPoints()
{
  for(unsigned int i = 0; i < Coordinates.size(); ++i)
    {
    this->CurrentRenderer->RemoveViewProp( Numbers[i]);
    this->CurrentRenderer->RemoveViewProp( Points[i]);
    }
  Numbers.clear();
  Points.clear();
  Coordinates.clear();
}
*/

void PointSelectionStyle2D::AddNumber(const double p[3])
{
  // Convert the current number to a string
  std::stringstream ss;
  ss << Coordinates.size();
  
  Coord3D coord;
  coord.x = p[0];
  coord.y = p[1];
  Coordinates.push_back(coord);

  // The coordinate provided is the corner of the pixel - we want to display the sphere in the center of the pixel.
  double markerCenter[3];
  markerCenter[0] = static_cast<int>( p[0] + 0.5 );
  markerCenter[1] = static_cast<int>( p[1] + 0.5 );
  markerCenter[2] = 0;
  std::cout << "Adding marker at " << markerCenter[0] << " "
            << markerCenter[1] << " " << markerCenter[2] << std::endl;

  // Create the number
  // Create the text
  vtkSmartPointer<vtkCaptionActor2D> captionActor = vtkSmartPointer<vtkCaptionActor2D>::New();
  captionActor->SetCaption( ss.str().c_str() );
  captionActor->SetAttachmentPoint(markerCenter);
  captionActor->BorderOff();
  //captionActor->SetPadding(10);
  captionActor->GetCaptionTextProperty()->SetFontFamilyToTimes();
  captionActor->GetCaptionTextProperty()->BoldOff();
  captionActor->GetCaptionTextProperty()->ItalicOff();
  captionActor->GetCaptionTextProperty()->ShadowOff();
  captionActor->ThreeDimensionalLeaderOff();
  this->Numbers.push_back(captionActor);
  this->CurrentRenderer->AddViewProp( captionActor );

  // Create the dot
  // Create a sphere
  vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
  sphereSource->SetRadius(.5);
  sphereSource->SetCenter(markerCenter);
  sphereSource->Update();

  // Create a mapper
  vtkSmartPointer<vtkPolyDataMapper> sphereMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  sphereMapper->SetInputConnection( sphereSource->GetOutputPort() );

  // Create a subclass of vtkActor: a vtkFollower that remains facing the camera
  vtkSmartPointer<vtkActor> sphereActor = vtkSmartPointer<vtkActor>::New();
  sphereActor->SetMapper( sphereMapper );
  sphereActor->GetProperty()->SetColor( 1, 0, 0 ); // red

  this->Points.push_back(sphereActor);
  //this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor( sphereActor );
  this->CurrentRenderer->AddViewProp( sphereActor );
}

void PointSelectionStyle2D::RemoveAll()
{
//   ResectioningHelpers::RemoveAllActors<vtkProp>(this->Numbers, this->CurrentRenderer);
//   ResectioningHelpers::RemoveAllActors<vtkActor>(this->Points, this->CurrentRenderer);
//   ResectioningHelpers::RemoveAllActors(this->Numbers, this->CurrentRenderer);
//   ResectioningHelpers::RemoveAllActors(this->Points, this->CurrentRenderer);

  this->Numbers.clear();
  this->Points.clear();
  this->Coordinates.clear();
}

void PointSelectionStyle2D::SetCurrentRenderer(vtkRenderer* const renderer)
{
  vtkInteractorStyleImage::SetCurrentRenderer(renderer);
}

void PointSelectionStyle2D::DeleteLastCorrespondence()
{
  this->CurrentRenderer->RemoveViewProp( this->Numbers[this->Numbers.size() - 1]);
  this->CurrentRenderer->RemoveViewProp( this->Points[this->Points.size() - 1]);
  this->Numbers.erase(this->Numbers.end()-1);
  this->Points.erase(this->Points.end()-1);
  this->Coordinates.erase(this->Coordinates.end()-1);
}
