#include "itkPTXReader.h"

#include <itkImage.h>
#include <itkImageFileWriter.h>

int main (int argc, char *argv[])
{
  if(argc != 3)
    {
    std::cout << "Required arguments: InputFilename(ptx) OutputFilename(mhd)" << std::endl;
    return EXIT_FAILURE;
    }

  std::string inputFilename = argv[1];
  std::string outputFilename = argv[2];

  itk::PTXReader::Pointer reader = itk::PTXReader::New();
  reader->SetFileName(inputFilename);
  reader->Update();

  reader->WriteRGBImage("test.png");

  typedef itk::ImageFileWriter< itk::Image<itk::CovariantVector<float, 8>, 2> > MetaWriterType;
  MetaWriterType::Pointer writer = MetaWriterType::New();
  writer->SetFileName(outputFilename);
  writer->SetInput(reader->GetOutput());
  writer->Update();

  return EXIT_SUCCESS;
}
