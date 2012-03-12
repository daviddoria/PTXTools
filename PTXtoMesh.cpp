#include "PTXImage.h"
#include "PTXReader.h"
#include "FilePrefix.h"

// VTK
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataWriter.h>

int main (int argc, char *argv[])
{
  if(argc != 3)
    {
    std::cout << "Required arguments: Input (PTX) OutputFileName" << std::endl;
    return EXIT_FAILURE;
    }

  std::string inputFilename = argv[1];
  std::string outputFileName = argv[2];

  //FilePrefix prefix(outputPrefix);
  

  PTXImage ptxImage = PTXReader::Read(inputFilename);

  vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
  ptxImage.GetMesh(polyData);

  vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName(outputFileName.c_str());
  writer->SetInputConnection(polyData->GetProducerPort());
  writer->Write();

  return EXIT_SUCCESS;
}
