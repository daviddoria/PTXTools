#include "Helpers.h"

#include <string>

#include <vtkVertexGlyphFilter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPolyData.h>

namespace Helpers
{

void OutputPolyData(vtkSmartPointer<vtkPolyData> points, std::string filename)
{
  // Output projected points for debugging
  vtkSmartPointer<vtkVertexGlyphFilter> vertexGlyphFilter =
    vtkSmartPointer<vtkVertexGlyphFilter>::New();
  vertexGlyphFilter->AddInput(points);
  vertexGlyphFilter->Update();

  vtkSmartPointer<vtkXMLPolyDataWriter> writer =
    vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName(filename.c_str());
  writer->SetInputConnection(vertexGlyphFilter->GetOutputPort());
  writer->Write();
}

}; // end Helpers namespace