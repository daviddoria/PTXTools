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

#include "Custom3DStyle.h"

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

vtkStandardNewMacro(Custom3DStyle);

Custom3DStyle::Custom3DStyle()
{

}
/*
void Custom3DStyle::Initialize()
{
  // This function must be called after calling selectionStyle->SetCurrentRenderer(renderer);
  if(!this->CurrentRenderer)
  {
    throw std::runtime_error("Initialize: CurrentRenderer is NULL!");
  }

  this->CurrentRenderer->AddViewProp(SelectedPointsActor);
  this->CurrentRenderer->AddViewProp(LabelActor);
}*/

void Custom3DStyle::OnLeftButtonDown()
{
  vtkPointPicker* picker = vtkPointPicker::SafeDownCast(this->Interactor->GetPicker());

  // If there is no picker, create one
  if(!picker)
  {
    picker = vtkPointPicker::New();
    Interactor->SetPicker(picker);
    std::cout << "Created picker." << std::endl;
  }

  int success = picker->Pick(this->Interactor->GetEventPosition()[0],
          this->Interactor->GetEventPosition()[1],
          0,  // always zero.
          this->CurrentRenderer);
  if(!success)
    {
    // Forward events
    vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
    return;
    }

  double picked[3] = {0,0,0};

  picker->GetPickPosition(picked);

  if(this->Interactor->GetShiftKey())
    {
    this->CurrentRenderer->GetActiveCamera()->SetFocalPoint(picked);
    }

  // Forward events
  vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
}
/*
void Custom3DStyle::SetCurrentRenderer(vtkRenderer* const renderer)
{
  vtkInteractorStyle::SetCurrentRenderer(renderer);
}*/
