#ifndef __itkPTXReader_h
#define __itkPTXReader_h

#include "itkImageSource.h"
#include "itkCovariantVector.h"

namespace itk
{
class PTXReader : public ImageSource<itk::Image<itk::CovariantVector< float, 8>, 2> >
{
  typedef itk::Image<itk::CovariantVector< float, 8>, 2> FullImageType;

public:
  /** Standard class typedefs. */
  typedef PTXReader             Self;
  typedef SmartPointer< Self >        Pointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(PTXReader, ImageSource);

  itkSetStringMacro(FileName);
  itkGetStringMacro(FileName);

  void WriteRGBImage(std::string filename);

  FullImageType::Pointer GetFullImage();

protected:
  PTXReader(){}
  ~PTXReader(){}

  /** Does the real work. */
  virtual void GenerateData();

  std::string m_FileName; // The file to be read


  FullImageType::Pointer FullImage;

private:
  PTXReader(const Self &); //purposely not implemented
  void operator=(const Self &);  //purposely not implemented

};
} //namespace ITK


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkPTXReader.txx"
#endif


#endif // __itkPTXReader_h