#ifndef ResectioningHelpers_HPP
#define ResectioningHelpers_HPP

namespace ResectioningHelpers
{
  
template<typename T>
//void RemoveAllActors(const std::vector<T>& actors, vtkRenderer* const renderer)
void RemoveAllActors(std::vector<T> actors, vtkRenderer* const renderer)
{
  for(unsigned int i = 0; i < actors.size(); ++i)
    {
    renderer->RemoveViewProp( actors[i]);
    }
}

template<typename TImage>
void DeepCopyScalarImage(const TImage* const input, TImage* const output)
{
  output->SetRegions(input->GetLargestPossibleRegion());
  output->Allocate();

  itk::ImageRegionConstIterator<TImage> inputIterator(input, input->GetLargestPossibleRegion());
  itk::ImageRegionIterator<TImage> outputIterator(output, output->GetLargestPossibleRegion());

  while(!inputIterator.IsAtEnd())
    {
    outputIterator.Set(inputIterator.Get());
    ++inputIterator;
    ++outputIterator;
    }
}

template<typename TImage>
void DeepCopyVectorImage(const TImage* const input, TImage* const output)
{
  output->SetRegions(input->GetLargestPossibleRegion());
  output->SetNumberOfComponentsPerPixel(input->GetNumberOfComponentsPerPixel());
  output->Allocate();

  itk::ImageRegionConstIterator<TImage> inputIterator(input, input->GetLargestPossibleRegion());
  itk::ImageRegionIterator<TImage> outputIterator(output, output->GetLargestPossibleRegion());

  while(!inputIterator.IsAtEnd())
    {
    outputIterator.Set(inputIterator.Get());
    ++inputIterator;
    ++outputIterator;
    }
}

} // end namespace

#endif
