#include "itkPTXRGBDIReader.h"

#include <itkImage.h>
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>

int main (int argc, char *argv[])
{
  if(argc != 3)
    {
    std::cout << "Required arguments: InputFilename(ptx) OutputFilename(mhd)" << std::endl;
    return EXIT_FAILURE;
    }

  std::string inputFilename = argv[1];
  std::string outputFilename = argv[2];

  itk::PTXRGBDIReader::Pointer reader = itk::PTXRGBDIReader::New();
  reader->SetFileName(inputFilename);
  reader->Update();

  typedef  itk::ImageFileWriter< itk::Image<itk::CovariantVector<float, 5>, 2> > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(outputFilename);
  writer->SetInput(reader->GetOutput());
  writer->Update();

  typedef  itk::ImageFileReader< itk::Image<itk::CovariantVector<float, 5>, 2> > ReaderType;
  ReaderType::Pointer metaReader = ReaderType::New();
  metaReader->SetFileName(outputFilename);
  metaReader->Update();

  return EXIT_SUCCESS;
}
