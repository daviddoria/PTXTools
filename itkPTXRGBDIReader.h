#ifndef __itkPTXRGBDIReader_h
#define __itkPTXRGBDIReader_h

#include "itkImageSource.h"
#include "itkCovariantVector.h"

namespace itk
{
class PTXRGBDIReader : public ImageSource<itk::Image<itk::CovariantVector< float, 5>, 2> >
{
public:
  /** Standard class typedefs. */
  typedef PTXRGBDIReader             Self;
  typedef SmartPointer< Self >        Pointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(PTXRGBDIReader, ImageSource);

  itkSetStringMacro(FileName);
  itkGetStringMacro(FileName);

protected:
  PTXRGBDIReader(){}
  ~PTXRGBDIReader(){}

  /** Does the real work. */
  virtual void GenerateData();

  std::string m_FileName; // The file to be read

private:
  PTXRGBDIReader(const Self &); //purposely not implemented
  void operator=(const Self &);  //purposely not implemented

};
} //namespace ITK


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkPTXRGBDIReader.txx"
#endif


#endif // __itkPTXRGBDIReader_h