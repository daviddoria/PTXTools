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
#include <vtkCamera.h>
#include <vtkFollower.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPointPicker.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRendererCollection.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkVectorText.h>

// STL
#include <sstream>
#include <stdexcept>

// Custom
#include "Helpers/Helpers.h"

vtkStandardNewMacro(PointSelectionStyle3D);

PointSelectionStyle3D::PointSelectionStyle3D()
{

}

void PointSelectionStyle3D::OnLeftButtonDown() 
{
  //if(!GetCurrentRenderer())
  if(!CurrentRenderer)
  {
    throw std::runtime_error("Must first call SetCurrentRenderer()!");
  }

  if(!Interactor)
  {
    throw std::runtime_error("Interactor is NULL!");
  }

  vtkPointPicker* picker = vtkPointPicker::SafeDownCast(this->Interactor->GetPicker());
  if(!picker)
  {
    // throw std::runtime_error("picker is NULL!");

    // If there is no picker, create one
    PointPicker = vtkSmartPointer<vtkPointPicker>::New();
    // PointPicker->PickFromListOn(); // If this is set to ON, then you must provide actors via AddPickList
    //PointPicker->AddPickList(pane->PointCloudActor);
    Interactor->SetPicker(PointPicker);

    picker = PointPicker;
  }

  //std::cout << "Picking pixel: " << this->Interactor->GetEventPosition()[0] << " " << this->Interactor->GetEventPosition()[1] << std::endl;
  int success = picker->Pick(this->Interactor->GetEventPosition()[0],
                             this->Interactor->GetEventPosition()[1],
                             0,  // always zero.
                             this->CurrentRenderer);
  if(!success)
    {
    // Forward events
    //std::cerr << "bad pick." << std::endl;
    vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
    return;
    }

  if(picker->GetDataSet() != this->Data)
    {
    //std::cerr << "Picked wrong data set: pointPicker->GetDataSet(): " << picker->GetDataSet() << " vs this->Data: " << this->Data << std::endl;
    vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
    return;
    //throw std::runtime_error("Did not pick from the correct data set!");
    }
  /*
  vtkIdType pointId = vtkPointPicker::SafeDownCast(this->Interactor->GetPicker())->GetPointId();
  double p[3];
  this->Data->GetPoint(pointId, p);
  
  //std::cout << "Picked point: " << pointId << std::endl;
  //std::cout << "Point: " << pointId << " should have coordinate: " << p[0] << " " << p[1] << " " << p[2] << std::endl;
  */
  double picked[3] = {0,0,0};
  
  picker->GetPickPosition(picked);
  //std::cout << "Picked point with coordinate: " << picked[0] << " " << picked[1] << " " << picked[2] << std::endl;

  if(this->Interactor->GetShiftKey())
    {
    this->CurrentRenderer->GetActiveCamera()->SetFocalPoint(picked);
    }

  // Only select the point if control is held
  if(this->Interactor->GetControlKey())
    {
    std::cout << "Picked 3D position " << picked[0] << " " << picked[1] << " " << picked[2] << std::endl;
    vtkIntArray* originalPixelArray = vtkIntArray::SafeDownCast(Data->GetPointData()->GetArray("OriginalPixel"));
    if(originalPixelArray)
      {
      int pixel[2];
      originalPixelArray->GetTupleValue(picker->GetPointId(), pixel);
      std::cout << "Picked PTX pixel " << pixel[0] << " " << pixel[1] << std::endl;
      }
    }

  // Forward events
  vtkInteractorStyleTrackballCamera::OnLeftButtonDown();

}

void PointSelectionStyle3D::SetData(vtkPolyData* const data)
{
  this->Data = data;
}
