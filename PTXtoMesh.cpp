#include "PTXImage.h"
#include "PTXReader.h"
#include "FilePrefix.h"

// VTK
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataWriter.h>

int main (int argc, char *argv[])
{
  if(argc != 4)
    {
    std::cout << "Required arguments: input.ptx MaxEdgeLength OutputFileName.vtp" << std::endl;
    return EXIT_FAILURE;
    }

  std::string inputFileName = argv[1];
  std::stringstream ss;
  ss << argv[2];
  float maxEdgeLength = 0.0f;
  ss >> maxEdgeLength;
  std::string outputFileName = argv[3];

  std::cout << "Input file: " << inputFileName << std::endl;
  std::cout << "maxEdgeLength: " << maxEdgeLength << std::endl;
  std::cout << "Output file: " << outputFileName << std::endl;
  

  PTXImage ptxImage = PTXReader::Read(inputFileName);
  ptxImage.ComputeMesh(maxEdgeLength);

  vtkPolyData* polyData = ptxImage.GetMesh();

  vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName(outputFileName.c_str());
  writer->SetInputData(polyData);
  writer->Write();

  return EXIT_SUCCESS;
}
