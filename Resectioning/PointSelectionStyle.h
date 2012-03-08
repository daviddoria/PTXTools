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

#ifndef PointSelectionStyle_H
#define PointSelectionStyle_H

#include <vtkInteractorStyle.h>

// STL
#include <vector>

// VTK
#include <vtkActor.h>
#include <vtkRenderer.h>

// Custom
#include "Coord.h"

// Define interaction style
class PointSelectionStyle
{
  public:
    virtual void OnLeftButtonDown() = 0;

    virtual void AddNumber(double p[3]) = 0;

    virtual void RemoveAll() = 0;

    virtual void DeleteLastCorrespondence() = 0;

    virtual void SetCurrentRenderer(vtkRenderer*) = 0;

    void SetMarkerRadius(const float radius);

    unsigned int GetNumberOfCorrespondences();

    Coord3D GetCorrespondence(const unsigned int correspondenceId);

    virtual void Initialize() = 0;
    
  protected:
    std::vector<vtkProp*> Numbers;
    std::vector<vtkActor*> Points;
    std::vector<Coord3D> Coordinates; // This should really be a more abstract Coord (with unknown dimension at this time)

    float MarkerRadius;
};

#endif
