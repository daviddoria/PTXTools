#include "PointSelectionStyle.h"

unsigned int PointSelectionStyle::GetNumberOfCorrespondences()
{
  return this->Numbers.size();
}

Coord3D PointSelectionStyle::GetCorrespondence(const unsigned int correspondenceId)
{
  return this->Coordinates[correspondenceId];
}

void PointSelectionStyle::SetMarkerRadius(const float radius)
{
  this->MarkerRadius = radius;
  //this->DotSource->SetRadius(this->MarkerRadius);
  //this->DotSource->Update();
}

