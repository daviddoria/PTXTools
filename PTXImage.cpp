#include "PTXImage.h"
#include "Frame.h"

// Submodules
#include <Helpers/Helpers.h>
#include <ITKHelpers/ITKHelpers.h>
#include <VTKHelpers/VTKHelpers.h>
#include <Mask.h>

// ITK
#include "itkAzimuthElevationToCartesianTransform.h"
#include "itkComposeImageFilter.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkConstantBoundaryCondition.h"
#include "itkCovariantVector.h"
#include "itkDerivativeImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkLaplacianRecursiveGaussianImageFilter.h"
#include "itkPoint.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkShrinkImageFilter.h"
#include "itkTileImageFilter.h"

// VTK
#include <vtkAppendPolyData.h>
#include <vtkCellArray.h>
#include <vtkDelaunay2D.h>
#include <vtkFloatArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPlaneSource.h>
#include <vtkLineSource.h>
#include <vtkMath.h>
#include <vtkPlane.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include <vtkTriangle.h>
#include <vtkUnsignedCharArray.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLStructuredGridWriter.h>

// STL
#include <fstream>
#include <sstream>
#include <string>

PTXImage::PTXImage() : Debug(false)
{
  // Create the main image
  this->FullImage = FullImageType::New();
  this->OriginalFullImage = FullImageType::New();

  this->Mesh = vtkSmartPointer<vtkPolyData>::New();
}

PTXImage::PTXImage(FullImageType* const fullImage) : PTXImage()
{
  this->FullImage = FullImageType::New();
  ITKHelpers::DeepCopy(fullImage, FullImage.GetPointer());
}

void PTXImage::SetDebug(const bool value)
{
  this->Debug = value;
}

PTXPixel PTXImage::GetPTXPixel(const itk::Index<2>& pixel) const
{

  if(this->FullImage->GetLargestPossibleRegion().IsInside(pixel))
    {
    return this->FullImage->GetPixel(pixel);
    }
  else
    {
    std::stringstream ss;
    ss << "Pixel " << pixel << " is not inside FullImage: " << this->FullImage->GetLargestPossibleRegion() << std::endl;
    throw std::runtime_error(ss.str());
    // Should never get here, this is just to avoid no return value warning.
    PTXPixel zeroPixel;
    //zeroPixel.Fill(2);
    return zeroPixel;
    }
}

void PTXImage::AppendPTXRight(const PTXImage& ptxImage)
{
  if(ptxImage.GetHeight() != this->GetHeight())
    {
    std::cerr << "Images to append side-by-side must be the same height!" << std::endl
              << "Current image is height " << this->GetHeight() << " while image to append is "
              << " height " << ptxImage.GetHeight() << std::endl;
    return;
    }
  // Tile the images side-by-side
  typedef itk::TileImageFilter< FullImageType, FullImageType > TileFilterType;

  TileFilterType::Pointer tileImageFilter = TileFilterType::New();

  // The following means "append to the right"
  itk::FixedArray< unsigned int, 2 > layout;
  layout[0] = 2;
  layout[1] = 0;
  tileImageFilter->SetLayout( layout );

  PTXPixel fillerValue;
  tileImageFilter->SetDefaultPixelValue( fillerValue );

  tileImageFilter->SetInput(0, this->FullImage);
  tileImageFilter->SetInput(1, ptxImage.GetFullImage());
  tileImageFilter->Update();

  ITKHelpers::DeepCopy(tileImageFilter->GetOutput(), this->FullImage.GetPointer());
}

unsigned int PTXImage::CountValidPoints() const
{
  unsigned int totalPoints = this->FullImage->GetLargestPossibleRegion().GetSize()[0] *
                                this->FullImage->GetLargestPossibleRegion().GetSize()[1];
  return  totalPoints - CountInvalidPoints();
}

unsigned int PTXImage::CountZeroPoints() const
{
  itk::ImageRegionConstIterator<FullImageType> fullImageIterator(this->FullImage,
                                                                 this->FullImage->GetLargestPossibleRegion());

  unsigned int numberOfZeroPoints = 0;
  while(!fullImageIterator.IsAtEnd())
    {
    PTXPixel fullPixel = fullImageIterator.Get();

    if(fullPixel.IsZero())
      {
      numberOfZeroPoints++;
      }

    ++fullImageIterator;
    }

  //std::cout << "There are " << numberOfInvalidPoints << " invalid points." << std::endl;
  return numberOfZeroPoints;
}

unsigned int PTXImage::CountInvalidPoints() const
{
  itk::ImageRegionConstIterator<FullImageType> fullImageIterator(this->FullImage,
                                                                 this->FullImage->GetLargestPossibleRegion());

  unsigned int numberOfInvalidPoints = 0;
  while(!fullImageIterator.IsAtEnd())
    {
    PTXPixel fullPixel = fullImageIterator.Get();

    if(!fullPixel.IsValid())
      {
      numberOfInvalidPoints++;
      }

    ++fullImageIterator;
    }

  //std::cout << "There are " << numberOfInvalidPoints << " invalid points." << std::endl;
  return numberOfInvalidPoints;
}

void PTXImage::WriteXYZ(const std::string& filePrefix) const
{
  typedef itk::ImageFileWriter< XYZImageType > XYZWriterType;
  XYZWriterType::Pointer xyzWriter = XYZWriterType::New();
  std::stringstream ssXYZ;
  ssXYZ << filePrefix << "_xyz.mha";
  xyzWriter->SetFileName(ssXYZ.str());
  xyzWriter->SetInput(GetXYZImage());
  xyzWriter->Update();
}

PTXImage::XYZImageType::Pointer PTXImage::GetXYZImage() const
{
  typedef itk::ComposeImageFilter<FloatImageType,
                              XYZImageType> ComposeCovariantVectorImageFilterType;

  ComposeCovariantVectorImageFilterType::Pointer composeFilter = ComposeCovariantVectorImageFilterType::New();

  composeFilter->SetInput1(GetCoordinateImage(0));
  composeFilter->SetInput2(GetCoordinateImage(1));
  composeFilter->SetInput3(GetCoordinateImage(2));
  composeFilter->Update();

  return composeFilter->GetOutput();
}

void PTXImage::WriteXYZLaplacian(const std::string& filePrefix) const
{
  typedef itk::ImageFileWriter< XYZImageType > XYZWriterType;
  XYZWriterType::Pointer xyzWriter = XYZWriterType::New();
  std::stringstream ssXYZ;
  ssXYZ << filePrefix << "_xyzLaplacian.mha";
  xyzWriter->SetFileName(ssXYZ.str());
  xyzWriter->SetInput(GetXYZLaplacian());
  xyzWriter->Update();
}

PTXImage::XYZImageType::Pointer PTXImage::GetXYZLaplacian() const
{
  typedef itk::ComposeImageFilter<FloatImageType,
                              XYZImageType> ComposeImageFilterType;

  ComposeImageFilterType::Pointer composeFilter = ComposeImageFilterType::New();

  composeFilter->SetInput1(GetLaplacian(0));
  composeFilter->SetInput2(GetLaplacian(1));
  composeFilter->SetInput3(GetLaplacian(2));
  composeFilter->Update();

  return composeFilter->GetOutput();
}

PTXImage::FloatImageType::Pointer PTXImage::GetCoordinateImage(const unsigned int dimension) const
{
  itk::ImageRegionConstIterator<FullImageType> imageIterator(this->FullImage, this->FullImage->GetLargestPossibleRegion());

  FloatImageType::Pointer image = FloatImageType::New();
  image->SetRegions(this->FullImage->GetLargestPossibleRegion());
  image->Allocate();
  image->FillBuffer(0);

  while(!imageIterator.IsAtEnd())
    {
    PTXPixel pixel = imageIterator.Get();
    image->SetPixel(imageIterator.GetIndex(), pixel.GetCoordinate(dimension));
    ++imageIterator;
    }

  return image;
}

void PTXImage::WriteFloatImage(const FloatImageType* const image, const std::string& filename) const
{
  typedef itk::ImageFileWriter< FloatImageType > FloatWriterType;
  FloatWriterType::Pointer writer = FloatWriterType::New();
  writer->SetFileName(filename);
  writer->SetInput(image);
  writer->Update();
}

void PTXImage::WriteX(const std::string& filePrefix) const
{
  std::stringstream ssX;
  ssX << filePrefix << "_x.mha";
  WriteFloatImage(GetXImage(), ssX.str());
}

PTXImage::FloatImageType::Pointer PTXImage::GetXImage() const
{
  return GetCoordinateImage(0);
}

void PTXImage::WriteY(const std::string& filePrefix) const
{
  std::stringstream ss;
  ss << filePrefix << "_y.mha";
  WriteFloatImage(GetYImage(), ss.str());
}

PTXImage::FloatImageType::Pointer PTXImage::GetYImage() const
{
  return GetCoordinateImage(1);
}

void PTXImage::WriteZ(const std::string& filePrefix) const
{
  std::stringstream ss;
  ss << filePrefix << "_z.mha";
  WriteFloatImage(GetZImage(), ss.str());
}

PTXImage::FloatImageType::Pointer PTXImage::GetZImage() const
{
  return GetCoordinateImage(2);
}

PTXImage::FloatImageType::Pointer PTXImage::GetLaplacian(const unsigned int dimension) const
{
#if 0
  typedef itk::DerivativeImageFilter<FloatImageType, FloatImageType >  DerivativeFilterType;

  // Create and setup a derivative filter
  DerivativeFilterType::Pointer derivativeFilter = DerivativeFilterType::New();
  derivativeFilter->SetInput( GetCoordinateImage(dimension) );
  derivativeFilter->SetDirection(dimension);
  derivativeFilter->Update();

  return derivativeFilter->GetOutput();

#endif

  typedef itk::LaplacianRecursiveGaussianImageFilter<
    FloatImageType, FloatImageType >  LaplacianFilterType;

  LaplacianFilterType::Pointer laplacianFilter = LaplacianFilterType::New();
  laplacianFilter->SetInput( GetCoordinateImage(dimension) );
  laplacianFilter->Update();

  return laplacianFilter->GetOutput();
}

void PTXImage::WriteXLaplacian(const std::string& filePrefix) const
{
  std::stringstream ss;
  ss << filePrefix << "_xLaplacian.mha";
  WriteFloatImage(GetXLaplacian(), ss.str());
}

PTXImage::FloatImageType::Pointer PTXImage::GetXLaplacian() const
{
  return GetLaplacian(0);
}

void PTXImage::WriteYLaplacian(const std::string& filePrefix) const
{
  std::stringstream ss;
  ss << filePrefix << "_yLaplacian.mha";
  WriteFloatImage(GetYLaplacian(), ss.str());
}

PTXImage::FloatImageType::Pointer PTXImage::GetYLaplacian() const
{
  return GetLaplacian(1);
}

void PTXImage::WriteZLaplacian(const std::string& filePrefix) const
{
  std::stringstream ss;
  ss << filePrefix << "_zLaplacian.mha";
  WriteFloatImage(GetZLaplacian(), ss.str());
}

PTXImage::FloatImageType::Pointer PTXImage::GetZLaplacian() const
{
  return GetLaplacian(2);
}

void PTXImage::WriteDepthThresholdMask(const std::string& filename, const float depthThreshold) const
{
  MaskImageType::Pointer mask = MaskImageType::New();
  mask->SetRegions(this->FullImage->GetLargestPossibleRegion());
  mask->Allocate();
  mask->FillBuffer(0);

  // Setup iterators
  itk::ImageRegionIterator<MaskImageType> maskIterator(mask, mask->GetLargestPossibleRegion());

  itk::ImageRegionConstIterator<FullImageType> fullImageIterator(this->FullImage, this->FullImage->GetLargestPossibleRegion());

  while(!fullImageIterator.IsAtEnd())
    {
    float depth = fullImageIterator.Get().GetDepth();

    if(depth > depthThreshold)
      {
      maskIterator.Set(0); // zero (black) means "valid"
      }
    else
      {
      maskIterator.Set(255); // nonzero means "unknown/invalid"
      }

    ++fullImageIterator;
    ++maskIterator;
    }

  std::stringstream ss;
  ss << filename << "_DepthThresholdMask.png";

  typedef  itk::ImageFileWriter<MaskImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(ss.str());
  writer->SetInput(mask);
  writer->Update();
}

void PTXImage::CreateValidityImage(MaskImageType* const image) const
{
  // Non-zero pixels in this image indicate valid pixels.
  image->SetRegions(this->FullImage->GetLargestPossibleRegion());
  image->Allocate();
  image->FillBuffer(0);

  // Setup iterators
  itk::ImageRegionIterator<MaskImageType> validityIterator(image, image->GetLargestPossibleRegion());

  itk::ImageRegionConstIterator<FullImageType> fullImageIterator(this->FullImage, this->FullImage->GetLargestPossibleRegion());

  while(!fullImageIterator.IsAtEnd())
    {
    PTXPixel fullPixel = fullImageIterator.Get();

    if(fullPixel.IsValid())
      {
      validityIterator.Set(255); // non-zero means "valid"
      }
    else
      {
      validityIterator.Set(0); // zero means "unknown/invalid"
      }

    ++fullImageIterator;
    ++validityIterator;
    }

}

void PTXImage::WriteInvalidMask(const std::string& filename) const
{
  MaskImageType::Pointer mask = MaskImageType::New();
  mask->SetRegions(this->FullImage->GetLargestPossibleRegion());
  mask->Allocate();
  mask->FillBuffer(0);

  // Setup iterators
  itk::ImageRegionIterator<MaskImageType> maskIterator(mask, mask->GetLargestPossibleRegion());

  itk::ImageRegionConstIterator<FullImageType> fullImageIterator(this->FullImage, this->FullImage->GetLargestPossibleRegion());

  while(!fullImageIterator.IsAtEnd())
    {
    PTXPixel fullPixel = fullImageIterator.Get();

    if(fullPixel.IsValid())
      {
      maskIterator.Set(0); // zero (black) means "valid"
      }
    else
      {
      maskIterator.Set(255); // nonzero means "unknown/invalid"
      }

    ++fullImageIterator;
    ++maskIterator;
    }

  std::stringstream ss;
  ss << filename << "_ValidityMask.png";

  typedef  itk::ImageFileWriter<MaskImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(ss.str());
  writer->SetInput(mask);
  writer->Update();
}

void PTXImage::CreateRGBImage(RGBImageType* const image) const
{
  // Setup the image
  image->SetRegions(this->FullImage->GetLargestPossibleRegion());
  image->Allocate();
  image->FillBuffer(itk::NumericTraits<RGBImageType::PixelType>::Zero);

  // Setup iterators
  itk::ImageRegionIterator<RGBImageType> rgbImageIterator(image, image->GetLargestPossibleRegion());
  rgbImageIterator.GoToBegin();

  itk::ImageRegionIterator<FullImageType> fullImageIterator(this->FullImage,
                                                            this->FullImage->GetLargestPossibleRegion());
  fullImageIterator.GoToBegin();

  // Traverse the full image, extracting the color information and saving it in the RGB image
  while(!rgbImageIterator.IsAtEnd())
    {
    PTXPixel fullPixel = fullImageIterator.Get();

    //itk::CovariantVector<unsigned char, 3> rgbPixel;
    RGBImageType::PixelType rgbPixel;

    if(fullPixel.IsValid())
      {
      /*
      rgbPixel[0] = fullPixel.R;
      rgbPixel[1] = fullPixel.G;
      rgbPixel[2] = fullPixel.B;
      */
      rgbPixel.SetRed(fullPixel.R);
      rgbPixel.SetGreen(fullPixel.G);
      rgbPixel.SetBlue(fullPixel.B);
      }
    else
      {
      // Set the pixel to bright green
      rgbPixel.SetRed(0);
      //rgbPixel[1] = 255;
      rgbPixel.SetGreen(0);
      rgbPixel.SetBlue(0);
      }

    rgbImageIterator.Set(rgbPixel);

    ++rgbImageIterator;
    ++fullImageIterator;
    }

}

void PTXImage::WriteEverything(const FilePrefix& filePrefix) const
{
  WritePTX(filePrefix);

  WriteRGBImage(filePrefix);

  WriteDepthImage(filePrefix);

  WritePointCloud(filePrefix);

}

void PTXImage::WriteRGBImage(const FilePrefix& filePrefix) const
{
  // This a convenience function which simply calls CreateRGBImage and then writes the result to a file
  RGBImageType::Pointer image = RGBImageType::New();
  CreateRGBImage(image);

  std::stringstream ss;
  ss << filePrefix.prefix << "_RGB.png";

  typedef  itk::ImageFileWriter< RGBImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(ss.str());
  writer->SetInput(image);
  writer->Update();
}

void PTXImage::CreatePointCloud(vtkPolyData* const pointCloud) const
{
  // Create point and color arrays
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

  vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
  colors->SetNumberOfComponents(3);
  colors->SetName("Colors");

  vtkSmartPointer<vtkFloatArray> depthArray = vtkSmartPointer<vtkFloatArray>::New();
  depthArray->SetNumberOfComponents(1);
  depthArray->SetName("Depths");

  // This array is used if you do any point set processing - it allows you to map the points pack
  // to their position in the PTX grid.
  vtkSmartPointer<vtkIntArray> originalPixelArray = vtkSmartPointer<vtkIntArray>::New();
  originalPixelArray->SetNumberOfComponents(2);
  originalPixelArray->SetName("OriginalPixel");

  vtkSmartPointer<vtkFloatArray> intensities = vtkSmartPointer<vtkFloatArray>::New();
  intensities->SetNumberOfComponents(1);
  intensities->SetName("Intensity");

  // Iterate through the full image extracting the coordinate and color information and adding them to their respective arrays
  itk::ImageRegionConstIteratorWithIndex<FullImageType> imageIterator(this->FullImage, this->FullImage->GetLargestPossibleRegion());

  while(!imageIterator.IsAtEnd())
    {
    PTXPixel pixel = imageIterator.Get();
    unsigned char rgb[3];
    rgb[0] = pixel.R;
    rgb[1] = pixel.G;
    rgb[2] = pixel.B;

    if(pixel.IsValid())
      {
      colors->InsertNextTupleValue(rgb);
      points->InsertNextPoint(pixel.X, pixel.Y, pixel.Z);
      intensities->InsertNextValue(pixel.Intensity);
      depthArray->InsertNextValue(pixel.GetDepth());

      int originalPixel[2];
      originalPixel[0] = imageIterator.GetIndex()[0];
      originalPixel[1] = imageIterator.GetIndex()[1];

      originalPixelArray->InsertNextTupleValue(originalPixel);
      }

    ++imageIterator;
    }

  // Combine the coordinates and colors into a PolyData
  vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints(points);
  polydata->GetPointData()->SetScalars(colors);
  polydata->GetPointData()->AddArray(depthArray);
  polydata->GetPointData()->AddArray(originalPixelArray);
  polydata->GetPointData()->AddArray(intensities);

  // Create a vertex at each point
  vtkSmartPointer<vtkVertexGlyphFilter> vertexGlyphFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
  vertexGlyphFilter->SetInputData(polydata);
  vertexGlyphFilter->Update();

  // Add the ptx image size to the field data
  vtkSmartPointer<vtkIntArray> imageSizeArray = vtkSmartPointer<vtkIntArray>::New();
  int imageSize[2] = {static_cast<int>(this->GetSize()[0]), static_cast<int>(this->GetSize()[1])};
  imageSizeArray->SetNumberOfComponents(2);
  imageSizeArray->SetName("ImageSize");
  imageSizeArray->InsertNextTupleValue(imageSize);
  polydata->GetFieldData()->AddArray(imageSizeArray);

  pointCloud->ShallowCopy(vertexGlyphFilter->GetOutput());
}


void PTXImage::CreateStructuredGrid(vtkSmartPointer<vtkStructuredGrid> structuredGrid) const
{
  int dimensions[3] = {static_cast<int>(this->GetWidth()), static_cast<int>(this->GetHeight()), 1};
  structuredGrid->SetDimensions(dimensions);

  unsigned int totalPoints = this->GetWidth() * this->GetHeight();

  // Create point and color arrays
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  points->SetNumberOfPoints(totalPoints);

  vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
  colors->SetNumberOfComponents(3);
  colors->SetNumberOfTuples(totalPoints);
  colors->SetName("Colors");

  vtkSmartPointer<vtkFloatArray> depthArray = vtkSmartPointer<vtkFloatArray>::New();
  depthArray->SetNumberOfComponents(1);
  depthArray->SetNumberOfTuples(totalPoints);
  depthArray->SetName("Depths");

  vtkSmartPointer<vtkFloatArray> intensities = vtkSmartPointer<vtkFloatArray>::New();
  intensities->SetNumberOfComponents(1);
  intensities->SetNumberOfTuples(totalPoints);
  intensities->SetName("Intensity");

  // Iterate through the full image extracting the coordinate and color information and adding them to their respective arrays
  itk::ImageRegionConstIteratorWithIndex<FullImageType> imageIterator(this->FullImage, this->FullImage->GetLargestPossibleRegion());

  while(!imageIterator.IsAtEnd())
    {
    int queryPoint[3] = {imageIterator.GetIndex()[0], imageIterator.GetIndex()[1], 0};
    vtkIdType pointId = vtkStructuredData::ComputePointId(dimensions, queryPoint);

    PTXPixel pixel = imageIterator.Get();

    if(pixel.IsValid())
      {
      unsigned char rgb[3];
      rgb[0] = pixel.R;
      rgb[1] = pixel.G;
      rgb[2] = pixel.B;
      colors->SetTupleValue(pointId, rgb);

      points->SetPoint(pointId, pixel.X, pixel.Y, pixel.Z);

      intensities->SetValue(pointId, pixel.Intensity);

      depthArray->SetValue(pointId, pixel.GetDepth());
      }
    else
      {
      unsigned char rgb[3] = {0,0,0};
      colors->SetTupleValue(pointId, rgb);

      points->SetPoint(pointId, 0,0,0);

      intensities->SetValue(pointId, 0);

      depthArray->SetValue(pointId, 0);

      structuredGrid->BlankPoint(pointId);
      }
    ++imageIterator;
    }

  // Combine the coordinates and colors into a vtkStructuredGrid
  structuredGrid->SetPoints(points);
  structuredGrid->GetPointData()->SetScalars(colors);
  structuredGrid->GetPointData()->AddArray(depthArray);
  structuredGrid->GetPointData()->AddArray(intensities);

  // Create a vertex at each point
//   vtkSmartPointer<vtkVertexGlyphFilter> vertexGlyphFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
//   vertexGlyphFilter->AddInput(polydata);
//   vertexGlyphFilter->Update();
}

void PTXImage::WritePointCloud(const std::string& fileName) const
{
  vtkSmartPointer<vtkPolyData> pointCloud = vtkSmartPointer<vtkPolyData>::New();
  CreatePointCloud(pointCloud);

  vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName(fileName.c_str());
  writer->SetInputData(pointCloud);
  writer->Write();
}

void PTXImage::WriteStructuredGrid(const std::string& fileName) const
{
  vtkSmartPointer<vtkStructuredGrid> structuredGrid =
    vtkSmartPointer<vtkStructuredGrid>::New();
  CreateStructuredGrid(structuredGrid);

  vtkSmartPointer<vtkXMLStructuredGridWriter> writer = vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
  writer->SetFileName(fileName.c_str());
  writer->SetInputData(structuredGrid);
  writer->Write();
  
}

void PTXImage::WritePointCloud(const FilePrefix& filePrefix) const
{
  // This a convenience function which simply calls CreatePointCloud and then writes the result to a file
  std::stringstream ss;
  ss << filePrefix.prefix << ".vtp";

  WritePointCloud(ss.str());
}

void PTXImage::CreateDepthImage(FloatImageType::Pointer image) const
{
  // Setup the image
  image->SetRegions(this->FullImage->GetLargestPossibleRegion());
  image->Allocate();
  image->FillBuffer(0);

  // Setup iterators
  itk::ImageRegionConstIteratorWithIndex<FullImageType> imageIterator(this->FullImage, this->FullImage->GetLargestPossibleRegion());
  itk::ImageRegionIterator<FloatImageType> depthImageIterator(image, image->GetLargestPossibleRegion());

  // Traverse the full image extracting the depth information and saving it in the depth image
  while(!imageIterator.IsAtEnd())
    {
    depthImageIterator.Set(imageIterator.Get().GetDepth());

    ++imageIterator;
    ++depthImageIterator;
    }

}

void PTXImage::WriteDepthImage(const FilePrefix& filePrefix) const
{
  // This a convenience function which simply calls CreateDepthImage and then writes the result to png (scaled) and mha (unscaled) files.
  // 'filePrefix' is, for example "filename" which will be used to internally produce the filenames "filename.mha" and "filename.png"

  FloatImageType::Pointer image = FloatImageType::New();
  CreateDepthImage(image);

  std::stringstream ssMHA;
  ssMHA << filePrefix.prefix << "Depth.mha";
  typedef  itk::ImageFileWriter<FloatImageType> MHAWriterType;
  MHAWriterType::Pointer mhaWriter = MHAWriterType::New();
  mhaWriter->SetFileName(ssMHA.str());
  mhaWriter->SetInput(image);
  mhaWriter->Update();

  typedef itk::RescaleIntensityImageFilter<FloatImageType, UnsignedCharImageType > RescaleFilterType;
  RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
  rescaleFilter->SetInput(image);
  rescaleFilter->SetOutputMinimum(0);
  rescaleFilter->SetOutputMaximum(255);
  rescaleFilter->Update();

  std::stringstream ssPNG;
  ssPNG << filePrefix.prefix << "Depth.png";

  typedef  itk::ImageFileWriter< UnsignedCharImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(ssPNG.str());
  writer->SetInput(rescaleFilter->GetOutput());
  writer->Update();

}

void PTXImage::WriteDepthLaplacian(const FilePrefix& filePrefix) const
{
  FloatImageType::Pointer image = GetDepthLaplacian();

  std::stringstream ss;
  ss << filePrefix.prefix << "_DepthLaplacian.mha";

  typedef  itk::ImageFileWriter< FloatImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(ss.str());
  writer->SetInput(image);
  writer->Update();
}

PTXImage::FloatImageType::Pointer PTXImage::GetDepthLaplacian() const
{
  FloatImageType::Pointer depthImage = FloatImageType::New();
  CreateDepthImage(depthImage);

  typedef itk::LaplacianRecursiveGaussianImageFilter<
    FloatImageType, FloatImageType >  LaplacianFilterType;
  LaplacianFilterType::Pointer laplacianFilter = LaplacianFilterType::New();
  laplacianFilter->SetInput( depthImage );
  laplacianFilter->Update();

  return laplacianFilter->GetOutput();
}

void PTXImage::CreateIntensityImage(FloatImageType::Pointer image) const
{
  // Setup the image
  image->SetRegions(this->FullImage->GetLargestPossibleRegion());
  image->Allocate();
  image->FillBuffer(0);

  // Setup iterators
  itk::ImageRegionConstIteratorWithIndex<FullImageType> imageIterator(this->FullImage, this->FullImage->GetLargestPossibleRegion());
  itk::ImageRegionIterator<FloatImageType> intensityImageIterator(image, image->GetLargestPossibleRegion());

  // Traverse the full image extracting the depth information and saving it in the depth image
  while(!imageIterator.IsAtEnd())
    {
    PTXPixel pixel = imageIterator.Get();
    if(pixel.IsValid())
      {
      intensityImageIterator.Set(pixel.Intensity);
      }
    else
      {
      intensityImageIterator.Set(0);
      }

    ++imageIterator;
    ++intensityImageIterator;
    }

}

void PTXImage::WriteIntensityImage(const FilePrefix& filePrefix) const
{
  // This a convenience function which simply calls CreateIntensityImage and then writes the result to png (scaled) and mha (unscaled) files.
  // 'filePrefix' is, for example "filename" which will be used to internally produce the filenames "filename.mha" and "filename.png"

  FloatImageType::Pointer image = FloatImageType::New();
  CreateIntensityImage(image);

  std::stringstream ssMHA;
  ssMHA << filePrefix.prefix << "Intensity.mha";
  typedef  itk::ImageFileWriter< FloatImageType > MHAWriterType;
  MHAWriterType::Pointer mhaWriter = MHAWriterType::New();
  mhaWriter->SetFileName(ssMHA.str());
  mhaWriter->SetInput(image);
  mhaWriter->Update();

  typedef itk::RescaleIntensityImageFilter< FloatImageType, UnsignedCharImageType > RescaleFilterType;
  RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
  rescaleFilter->SetInput(image);
  rescaleFilter->SetOutputMinimum(0);
  rescaleFilter->SetOutputMaximum(255);
  rescaleFilter->Update();

  std::stringstream ssPNG;
  ssPNG << filePrefix.prefix << "Intensity.png";

  typedef  itk::ImageFileWriter< UnsignedCharImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(ssPNG.str());
  writer->SetInput(rescaleFilter->GetOutput());
  writer->Update();

}

void PTXImage::ReplaceDepth(const FloatImageType* const depthImage)
{
  // This function allows the depth map to be modified externally and the new map applied to the grid

  if(OriginalFullImage->GetLargestPossibleRegion() != FullImage->GetLargestPossibleRegion())
  {
    throw std::runtime_error("You must call Backup() on the originally read ptx before using ReplaceDepth()!");
  }

  PTXImage originalPTXImage(OriginalFullImage);
  originalPTXImage.ComputeAverageDeltaPhi();
  originalPTXImage.ComputeAverageDeltaTheta();

  std::cout << "Original info:" << std::endl;
  originalPTXImage.OutputInfo();
  
  // Setup iterators
  itk::ImageRegionIteratorWithIndex<FullImageType> imageIterator(this->FullImage,
                                                                 this->FullImage->GetLargestPossibleRegion());
  itk::ImageRegionConstIterator<FloatImageType> newDepthIterator(depthImage,
                                                                 depthImage->GetLargestPossibleRegion());

  itk::Point<float, 3> origin;
  origin.Fill(0);

  unsigned int numberOfApproximatedDirections = 0;
  unsigned int numberOfInvalidPoints = 0;
  while(!imageIterator.IsAtEnd())
    {
    // Get the old point
    PTXPixel currentPTXPixel = imageIterator.Get();

    PTXPixel originalPTXPixel = OriginalFullImage->GetPixel(imageIterator.GetIndex());

    // If the point is not valid, set it to the origin and move on to the next point.
    if(!currentPTXPixel.Valid) // Note: this should not be IsValid(), as we test the IsZero condition separately
    {
      numberOfInvalidPoints++;
      currentPTXPixel.X = origin[0];
      currentPTXPixel.Y = origin[1];
      currentPTXPixel.Z = origin[2];

      imageIterator.Set(currentPTXPixel);
      ++imageIterator;
      ++newDepthIterator;
      continue;
    }

    // If the original point was valid, use it's position as the unit vector along which to set the new depth.
    // If the original point was not valid, interpolate this point's direction.
    itk::Vector<float, 3> rayDirection;

    if(originalPTXPixel.IsValid()) // This means the point is marked as valid, and is not zero
    {
      itk::Point<float, 3> oldPoint;
      oldPoint[0] = originalPTXPixel.X;
      oldPoint[1] = originalPTXPixel.Y;
      oldPoint[2] = originalPTXPixel.Z;
      //std::cout << "Old point " << oldPoint << std::endl;
      // For testing only, compute the old depth
      double oldDepth = origin.EuclideanDistanceTo(oldPoint);
      if(this->Debug)
        {
        std::cout << "Old depth: " << oldDepth << " New depth: " << newDepthIterator.Get() << std::endl;
        }

      //std::cout << "Actual theta: " << GetTheta(imageIterator.GetIndex())
      //          << " phi: " << GetPhi(imageIterator.GetIndex()) << std::endl;
      
      // Get the vector from the origin (scanner location) to the old point
      rayDirection = oldPoint - origin;
    }
    else
    {
      // We should get here if the point is 0 0 0 1 0 0 0 - indicating that it should be filled,
      // but was not previously valid.
      // std::cout << "Approximating ray direction..." << std::endl;
      numberOfApproximatedDirections++;
      
      rayDirection = originalPTXImage.ApproximateRayDirection(imageIterator.GetIndex());

//       std::cout << "Approximated ray direction " << rayDirection << std::endl;
//       std::cout << "New depth: " << newDepthIterator.Get() << std::endl;
      
      if(isnan(rayDirection[0]))
      {
        throw std::runtime_error("ReplaceDepth(): Error computing ray direction!");
      }
    } // end if(pixel is valid)

    // Get a unit vector in the direction of the old point
    rayDirection.Normalize();

    //std::cout << "Using ray direction: " << rayDirection << std::endl;
    // Compute the new point from the vector and the new depth
    itk::Point<float, 3> newPoint = origin + rayDirection * newDepthIterator.Get();

    // Save the new point in the PTXPixel
    currentPTXPixel.X = newPoint[0];
    currentPTXPixel.Y = newPoint[1];
    currentPTXPixel.Z = newPoint[2];

    imageIterator.Set(currentPTXPixel);

    ++imageIterator;
    ++newDepthIterator;
    }

  std::cout << "ReplaceDepth(): Approximated " << numberOfApproximatedDirections << " directions." << std::endl;
  std::cout << "ReplaceDepth(): Invalid points " << numberOfInvalidPoints << std::endl;
  
}

void PTXImage::ReplaceRGB(const RGBVectorImageType* const rgb)
{
  // Recolor valid points to match the input 'rgb'

  if(rgb->GetLargestPossibleRegion() != this->FullImage->GetLargestPossibleRegion())
    {
    std::stringstream ss;
    ss << "PTXImage::ReplaceRGB: Input image must be exactly the same size as the PTX file!" << std::endl
       << "Input image is " << rgb->GetLargestPossibleRegion()
       << " and PTX is " << this->FullImage->GetLargestPossibleRegion() << std::endl;
    throw std::runtime_error(ss.str());
    }

  // Setup iterators
  itk::ImageRegionIterator<FullImageType> imageIterator(this->FullImage,
                                                        this->FullImage->GetLargestPossibleRegion());
  itk::ImageRegionConstIterator<RGBVectorImageType> rgbIterator(rgb, rgb->GetLargestPossibleRegion());

  while(!imageIterator.IsAtEnd())
    {
    // Get the old point
    PTXPixel pixel = imageIterator.Get();
    if(pixel.IsValid())
      {
      // Copy the color from the RGB image
      pixel.R = rgbIterator.Get()[0];
      pixel.G = rgbIterator.Get()[1];
      pixel.B = rgbIterator.Get()[2];
      imageIterator.Set(pixel);
      }

    ++imageIterator;
    ++rgbIterator;
    }
}


void PTXImage::ReplaceRGB(const RGBImageType* const rgb)
{
  // Recolor valid points to match the input 'rgb'

  if(rgb->GetLargestPossibleRegion() != this->FullImage->GetLargestPossibleRegion())
    {
    std::cerr << "Input image must be exactly the same size as the PTX file!" << std::endl;
    std::cerr << "Input image is " << rgb->GetLargestPossibleRegion()
              << " and PTX is " << this->FullImage->GetLargestPossibleRegion() << std::endl;
    return;
    }

  // Setup iterators
  itk::ImageRegionIterator<FullImageType> imageIterator(this->FullImage,
                                                        this->FullImage->GetLargestPossibleRegion());
  itk::ImageRegionConstIterator<RGBImageType> rgbIterator(rgb, rgb->GetLargestPossibleRegion());

  while(!imageIterator.IsAtEnd())
    {
    // Get the old point
    PTXPixel pixel = imageIterator.Get();
    if(pixel.IsValid())
      {
      // Copy the color from the RGB image
      pixel.R = rgbIterator.Get().GetRed();
      pixel.G = rgbIterator.Get().GetGreen();
      pixel.B = rgbIterator.Get().GetBlue();
      imageIterator.Set(pixel);
      }

    ++imageIterator;
    ++rgbIterator;
    }
}

void PTXImage::SetAllPointsToValid()
{
  itk::ImageRegionIteratorWithIndex<FullImageType> imageIterator(this->FullImage,
                                                                 this->FullImage->GetLargestPossibleRegion());

  while(!imageIterator.IsAtEnd())
    {
    PTXPixel pixel = imageIterator.Get();
    pixel.Valid = true;
    imageIterator.Set(pixel);
    ++imageIterator;
    }
}

void PTXImage::ReplaceValidity(const MaskImageType* const validityImage)
{
  if(validityImage->GetLargestPossibleRegion() != this->FullImage->GetLargestPossibleRegion())
    {
    std::cerr << "ReplaceValidity(): Input image must be exactly the same size as the PTX file!" << std::endl;
    std::cerr << "Input image is " << validityImage->GetLargestPossibleRegion()
              << " and PTX is " << this->FullImage->GetLargestPossibleRegion() << std::endl;
    return;
    }

  itk::ImageRegionIteratorWithIndex<FullImageType> imageIterator(this->FullImage,
                                                                 this->FullImage->GetLargestPossibleRegion());

  while(!imageIterator.IsAtEnd())
    {
    PTXPixel pixel = imageIterator.Get();
    pixel.Valid = validityImage->GetPixel(imageIterator.GetIndex());
    imageIterator.Set(pixel);
    ++imageIterator;
    }
}

// void PTXImage::ReplaceRGB(const RGBImageType::Pointer rgb)
// {
//   // Setup iterators
//   itk::ImageRegionIterator<FullImageType> fullImageIterator(this->FullImage,
//                                                             this->FullImage->GetLargestPossibleRegion());
//   itk::ImageRegionConstIterator<RGBImageType> rgbIterator(rgb, rgb->GetLargestPossibleRegion());
//
//   while(!fullImageIterator.IsAtEnd())
//     {
//     // Get the old point
//     PTXPixel pixel = fullImageIterator.Get();
//
//     // Copy the color from the RGB image
//     pixel.R = rgbIterator.Get().GetRed();
//     pixel.G = rgbIterator.Get().GetGreen();
//     pixel.B = rgbIterator.Get().GetBlue();
//     fullImageIterator.Set(pixel);
//
//     if(pixel.Valid)
//       {
//
//       }
//
//     ++fullImageIterator;
//     ++rgbIterator;
//     }
//
// }


void PTXImage::ReplaceXYZ(const XYZImageType* const xyz)
{
  // Setup iterators
  itk::ImageRegionIterator<FullImageType> imageIterator(this->FullImage,
                                                        this->FullImage->GetLargestPossibleRegion());
  itk::ImageRegionConstIterator<XYZImageType> xyzIterator(xyz, xyz->GetLargestPossibleRegion());

  unsigned int newPoints = 0; // These points were previously invalid. This is just record keeping for fun.

  while(!imageIterator.IsAtEnd())
    {
    // Get the old point
    PTXPixel pixel = imageIterator.Get();

    // Copy the point from the XYZ image
    pixel.X = xyzIterator.Get()[0];
    pixel.Y = xyzIterator.Get()[1];
    pixel.Z = xyzIterator.Get()[2];
    imageIterator.Set(pixel);

    if(pixel.Valid)
      {
      newPoints++;
      }

    ++imageIterator;
    ++xyzIterator;
    }

  std::cout << "There were " << newPoints << " new points (previously invalid)" << std::endl;

}

void PTXImage::ReplaceRGBD(const RGBDImageType* const rgbdImage)
{
  // Replace RGB
  std::vector<unsigned int> rgbChannels = {0, 1, 2};
  RGBVectorImageType::Pointer rgbImage = RGBVectorImageType::New();
  ITKHelpers::ExtractChannels(rgbdImage, rgbChannels, rgbImage.GetPointer());

  ReplaceRGB(rgbImage.GetPointer());

  // Replace Depth
  DepthImageType::Pointer depthImage = DepthImageType::New();
  ITKHelpers::ExtractChannel(rgbdImage, 3, depthImage.GetPointer());

  ReplaceDepth(depthImage.GetPointer());
}


itk::Vector<float, 3> PTXImage::ApproximateRayDirection(const itk::Index<2>& queryPixel) const
{
  /**
   * This function creates a point unit distance from the origin in the approximate direction
   * the old point would have been aquired. It should be used, for example, when replacing a
   * previously invalid point.
   */
  float phi = ApproximatePhi(queryPixel);
  float theta = ApproximateTheta(queryPixel);

  // std::cout << "Approximated Theta: " << theta << " phi: " << phi << std::endl;

//   if(isnan(phi))
//   {
//     throw std::runtime_error("Phi is nan!");
//   }
// 
//   if(isnan(theta))
//   {
//     //ApproximateTheta(pixel);
//     throw std::runtime_error("Theta is nan!");
//   }

  double x, y, z;
  Helpers::SphericalToCartesian(x,y,z,1.0f,theta,phi);

  itk::Vector<float, 3> direction;
  direction[0] = x;
  direction[1] = y;
  direction[2] = z;
  direction.Normalize();
//  std::cout << "Approximated direction: " << direction << std::endl;
  
  return direction;
}

itk::Point<float, 3> PTXImage::ApproximateOldPoint(const itk::Index<2>& pixel) const
{
  /**
   * This function creates a point unit distance from the origin in the approximate direction
   * the old point would have been aquired. It should be used, for example, when replacing a
   * previously invalid point.
   */
  float phi = ApproximatePhi(pixel);
  float theta = ApproximateTheta(pixel);

  if(isnan(phi))
  {
    throw std::runtime_error("Phi is nan!");
  }

  if(isnan(theta))
  {
    ApproximateTheta(pixel);
    throw std::runtime_error("Theta is nan!");
  }
  
  typedef itk::Point<float, 3> PointType;
  PointType azEl;
  azEl[0] = theta;
  azEl[1] = phi;
  azEl[2] = 1;

  typedef itk::AzimuthElevationToCartesianTransform< float, 3 >
    AzimuthElevationToCartesian;
  AzimuthElevationToCartesian::Pointer azimuthElevation =
    AzimuthElevationToCartesian::New();
  return azimuthElevation->TransformAzElToCartesian(azEl);
}

void PTXImage::ApplyMask(const MaskImageType* const mask)
{
  // Blank the PTX image in areas where mask is non-zero
  itk::ImageRegionIterator<FullImageType> imageIterator(this->FullImage, this->FullImage->GetLargestPossibleRegion());
  itk::ImageRegionConstIterator<MaskImageType> maskIterator(mask, mask->GetLargestPossibleRegion());

  while(!imageIterator.IsAtEnd())
    {
    if(maskIterator.Get())
      {
      PTXPixel pixel = imageIterator.Get();
      pixel.Valid = false;
      imageIterator.Set(pixel);
      }
    ++imageIterator;
    ++maskIterator;
    }
}

void PTXImage::CreateRGBDIImage(RGBDIImageType::Pointer image) const
{
  // Create the 5 channels
  FloatImageType::Pointer intensityImage = FloatImageType::New();
  CreateIntensityImage(intensityImage);

  FloatImageType::Pointer depthImage = FloatImageType::New();
  CreateDepthImage(depthImage);

  RGBImageType::Pointer rgbImage = RGBImageType::New();
  CreateRGBImage(rgbImage);

  // Setup the image
  image->SetRegions(this->FullImage->GetLargestPossibleRegion());
  image->Allocate();
  image->FillBuffer(itk::NumericTraits<RGBDIImageType::PixelType>::Zero);

  // Setup iterators
  itk::ImageRegionIterator<RGBDIImageType> imageIterator(image, image->GetLargestPossibleRegion());
  itk::ImageRegionConstIterator<FloatImageType> intensityImageIterator(intensityImage, intensityImage->GetLargestPossibleRegion());
  itk::ImageRegionConstIterator<FloatImageType> depthImageIterator(depthImage, depthImage->GetLargestPossibleRegion());
  itk::ImageRegionConstIterator<RGBImageType> rgbImageIterator(rgbImage, rgbImage->GetLargestPossibleRegion());

  // Copy the images into their respective channels
  while(!imageIterator.IsAtEnd())
    {

    RGBDIImageType::PixelType pixel;
    pixel[0] = rgbImageIterator.Get()[0];
    pixel[1] = rgbImageIterator.Get()[1];
    pixel[2] = rgbImageIterator.Get()[2];
    pixel[3] = depthImageIterator.Get();
    pixel[4] = intensityImageIterator.Get();

    imageIterator.Set(pixel);

    ++imageIterator;
    ++intensityImageIterator;
    ++depthImageIterator;
    ++rgbImageIterator;
    }
}

void PTXImage::WriteRGBDIImage(const FilePrefix& filePrefix) const
{
  RGBDIImageType::Pointer image = RGBDIImageType::New();
  CreateRGBDIImage(image);

  std::stringstream ss;
  ss << filePrefix.prefix << "_RGBDI.mha";

  typedef  itk::ImageFileWriter< RGBDIImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(ss.str());
  writer->SetInput(image);
  writer->Update();
}

void PTXImage::CreateRGBDImage(RGBDImageType::Pointer image) const
{
  // Create the 4 channels
  FloatImageType::Pointer depthImage = FloatImageType::New();
  CreateDepthImage(depthImage);

  RGBImageType::Pointer rgbImage = RGBImageType::New();
  CreateRGBImage(rgbImage);

  // Setup the image
  image->SetRegions(this->FullImage->GetLargestPossibleRegion());
  image->Allocate();
  image->FillBuffer(itk::NumericTraits<RGBDImageType::PixelType>::Zero);

  // Setup iterators
  itk::ImageRegionIterator<RGBDImageType> imageIterator(image, image->GetLargestPossibleRegion());
  itk::ImageRegionConstIterator<FloatImageType> depthImageIterator(depthImage, depthImage->GetLargestPossibleRegion());
  itk::ImageRegionConstIterator<RGBImageType> rgbImageIterator(rgbImage, rgbImage->GetLargestPossibleRegion());

  // Copy the images into their respective channels
  while(!imageIterator.IsAtEnd())
    {

    RGBDImageType::PixelType pixel;
    pixel[0] = rgbImageIterator.Get()[0];
    pixel[1] = rgbImageIterator.Get()[1];
    pixel[2] = rgbImageIterator.Get()[2];
    pixel[3] = depthImageIterator.Get();

    imageIterator.Set(pixel);

    ++imageIterator;
    ++depthImageIterator;
    ++rgbImageIterator;
    }
}

void PTXImage::WriteRGBDImage(const FilePrefix& filePrefix) const
{
  RGBDImageType::Pointer image = RGBDImageType::New();
  CreateRGBDImage(image);

  std::stringstream ss;
  ss << filePrefix.prefix << "_RGBD.mha";

  typedef  itk::ImageFileWriter< RGBDImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(ss.str());
  writer->SetInput(image);
  writer->Update();
}


void PTXImage::CreateRGBDVImage(RGBDVImageType::Pointer image) const
{
  // Create the 5 channels
  FloatImageType::Pointer depthImage = FloatImageType::New();
  CreateDepthImage(depthImage);

  RGBImageType::Pointer rgbImage = RGBImageType::New();
  CreateRGBImage(rgbImage);

  MaskImageType::Pointer validityImage = MaskImageType::New();
  CreateValidityImage(validityImage);

  // Setup the image
  image->SetRegions(this->FullImage->GetLargestPossibleRegion());
  image->Allocate();
  image->FillBuffer(itk::NumericTraits<RGBDVImageType::PixelType>::Zero);

  // Setup iterators
  itk::ImageRegionIterator<RGBDVImageType> imageIterator(image, image->GetLargestPossibleRegion());
  itk::ImageRegionConstIterator<FloatImageType> depthImageIterator(depthImage, depthImage->GetLargestPossibleRegion());
  itk::ImageRegionConstIterator<MaskImageType> validityImageIterator(validityImage, validityImage->GetLargestPossibleRegion());
  itk::ImageRegionConstIterator<RGBImageType> rgbImageIterator(rgbImage, rgbImage->GetLargestPossibleRegion());

  // Copy the images into their respective channels
  while(!imageIterator.IsAtEnd())
    {
    RGBDVImageType::PixelType pixel;
    pixel[0] = rgbImageIterator.Get()[0];
    pixel[1] = rgbImageIterator.Get()[1];
    pixel[2] = rgbImageIterator.Get()[2];
    pixel[3] = depthImageIterator.Get();
    pixel[4] = validityImageIterator.Get();

    imageIterator.Set(pixel);

    ++imageIterator;
    ++depthImageIterator;
    ++rgbImageIterator;
    ++validityImageIterator;
    }
}

void PTXImage::WriteRGBDVImage(const FilePrefix& filePrefix) const
{
  RGBDVImageType::Pointer image = RGBDVImageType::New();
  CreateRGBDVImage(image);

  std::stringstream ss;
  ss << filePrefix.prefix << "_RGBDV.mha";

  typedef  itk::ImageFileWriter< RGBDVImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(ss.str());
  writer->SetInput(image);
  writer->Update();
}

void PTXImage::ComputeAverageDeltaTheta()
{
  itk::ImageRegion<2> region = this->FullImage->GetLargestPossibleRegion();

  std::vector<float> deltaThetas;

  itk::Index<2> lastIndex;
  lastIndex[0] = -1;
  lastIndex[1] = -1;

  for(unsigned int row = 0; row < region.GetSize()[1]; row++)
    {
    for(unsigned int col = 0; col < region.GetSize()[0]; col++)
      {
      itk::Index<2> index;
      index[0] = col;
      index[1] = row;

      if(!region.IsInside(lastIndex))
      {
        lastIndex = index;
        continue;
      }

      if(!this->FullImage->GetPixel(index).IsValid() || // the current pixel is invalid
        !this->FullImage->GetPixel(lastIndex).IsValid() || // the previous pixel was invalid
        this->FullImage->GetPixel(index).IsZero() || // this pixel is marked as valid, but is (0,0,0)
        this->FullImage->GetPixel(lastIndex).IsZero() || // this pixel is marked as valid, but is (0,0,0)
        lastIndex[1] != index[1] // if we are not on the same row as the last pixel
        )
        {
        // do nothing
        }
      else
        {
        deltaThetas.push_back(GetTheta(index) - GetTheta(lastIndex));
        }

      lastIndex = index;
      }
    }

  float total = 0.0;
  for(unsigned int i = 0; i < deltaThetas.size(); i++)
    {
    total += fabs(deltaThetas[i]);
    }

  if(deltaThetas.size() == 0)
    {
    throw std::runtime_error("There are no valid theta neighbors!");
    }

  this->AverageDeltaTheta = total/static_cast<float>(deltaThetas.size());

  std::cout << "AverageDeltaTheta: " << this->AverageDeltaTheta << std::endl;
}

void PTXImage::ComputeAverageDeltaPhi()
{
  itk::ImageRegion<2> region = this->FullImage->GetLargestPossibleRegion();

  std::vector<float> deltaPhis;

  itk::Index<2> lastIndex;
  lastIndex[0] = -1;
  lastIndex[1] = -1;

  for(unsigned int col = 0; col < region.GetSize()[0]; col++)
    {
    for(unsigned int row = 0; row < region.GetSize()[1]; row++)
      {
      itk::Index<2> index;
      index[0] = col;
      index[1] = row;

      if(!region.IsInside(lastIndex))
      {
        lastIndex = index;
        continue;
      }
      
      if(!this->FullImage->GetPixel(index).IsValid() || // the current pixel is invalid
        !this->FullImage->GetPixel(lastIndex).IsValid() || // the previous pixel is invalid
        this->FullImage->GetPixel(index).IsZero() || // this pixel is marked as valid, but is (0,0,0)
        this->FullImage->GetPixel(lastIndex).IsZero() || // this pixel is marked as valid, but is (0,0,0)
        lastIndex[0] != index[0] // if we are not on the same column as the last pixel
        )
        {
        // do nothing
        }
      else
        {
        deltaPhis.push_back(fabs(GetPhi(index) - GetPhi(lastIndex)));
        }

      lastIndex = index;
      }
    }

  float total = 0.0;
  for(unsigned int i = 0; i < deltaPhis.size(); i++)
    {
    total += fabs(deltaPhis[i]);
    }

  if(deltaPhis.size() == 0)
    {
    throw std::runtime_error("There are no valid phi neighbors!");
    }

  this->AverageDeltaPhi = total/static_cast<float>(deltaPhis.size());

  std::cout << "AverageDeltaPhi: " << this->AverageDeltaPhi << std::endl;
}

float PTXImage::GetPhi(const itk::Index<2>& index) const
{
  if(this->FullImage->GetPixel(index).IsZero())
  {
    throw std::runtime_error("GetPhi index has a zero point!");
  }
  float x = this->FullImage->GetPixel(index).X;
  float y = this->FullImage->GetPixel(index).Y;
  float z = this->FullImage->GetPixel(index).Z;
  
  double r, theta, phi;
  
  Helpers::CartesianToSpherical(r, theta, phi, x, y, z );

  return phi;
}

float PTXImage::GetTheta(const itk::Index<2>& index) const
{
  if(this->FullImage->GetPixel(index).IsZero())
  {
    throw std::runtime_error("GetTheta index has a zero point!");
  }
  
  float x = this->FullImage->GetPixel(index).X;
  float y = this->FullImage->GetPixel(index).Y;
  float z = this->FullImage->GetPixel(index).Z;

  double r, theta, phi;

  Helpers::CartesianToSpherical(r, theta, phi, x, y, z );

  return theta;
}

float PTXImage::ApproximateTheta(const itk::Index<2>& queryPixel) const
{
  itk::Index<2> corner1 = {{0,0}};
//  PrintCoordinate(corner1);
  
  itk::Index<2> corner2 = {{static_cast<itk::Index<2>::IndexValueType>(this->GetFullRegion().GetSize()[0] - 1), 0}};
//  PrintCoordinate(corner2);

  itk::Index<2> left = corner2;
  itk::Index<2> right = corner1;
  
  if(this->FullImage->GetPixel(left).IsZero() || this->FullImage->GetPixel(right).IsZero())
  {
    throw std::runtime_error("ApproximateTheta: one of the corner is a zero point!");
  }

  if(this->Debug)
  {
    std::cout << "Left point: ";
    PrintCoordinate(left);

    std::cout << "Right point: ";
    PrintCoordinate(right);
  }

  float minTheta = GetTheta(left);
  float maxTheta = GetTheta(right);

  float thetaStep = fabs(maxTheta - minTheta)/ static_cast<float>(this->GetFullRegion().GetSize()[0]);
//  std::cout << "minTheta: " << minTheta << " maxTheta " << maxTheta << " thetaStep: " << thetaStep << std::endl;
  //float approximateTheta = minTheta + static_cast<float>(queryPixel[0]) * thetaStep;
  float approximateTheta = minTheta + static_cast<float>(this->GetFullRegion().GetSize()[0] - queryPixel[0]) * thetaStep;
//  std::cout << "approximateTheta: " << approximateTheta << std::endl;
  return approximateTheta;
}

float PTXImage::ApproximatePhi(const itk::Index<2>& queryPixel) const
{
//   itk::Index<2> bottom = {{0,0}};
//   itk::Index<2> top = {{0, maxIndex}};

  itk::Index<2> corner1 = {{0,0}};
//  PrintCoordinate(corner1);

  itk::Index<2> corner2 = {{0, static_cast<itk::Index<2>::IndexValueType>(this->GetFullRegion().GetSize()[1] - 1)}};
//  PrintCoordinate(corner2);
  
  itk::Index<2> bottom = corner2;
  itk::Index<2> top = corner1;

  if(this->FullImage->GetPixel(bottom).IsZero() || this->FullImage->GetPixel(top).IsZero())
  {
    throw std::runtime_error("ApproximateTheta: one of the corner is a zero point!");
  }
  
  float minPhi = GetPhi(bottom);
  float maxPhi = GetPhi(top);

  float phiStep = fabs(maxPhi - minPhi)/ static_cast<float>(this->GetFullRegion().GetSize()[1]);
  
//  std::cout << "minPhi: " << minPhi << " maxPhi " << maxPhi << " phiStep: " << phiStep << std::endl;
  
  float approximatePhi = minPhi + static_cast<float>(this->GetFullRegion().GetSize()[1] - queryPixel[1]) * phiStep;
  //float approximatePhi = minPhi + static_cast<float>(queryPixel[1]) * phiStep;
//  std::cout << "approximatePhi: " << approximatePhi << std::endl;
  
  return approximatePhi;
}

itk::Index<2> PTXImage::FindNearestValidPixel(const itk::Index<2>& pixel, itk::Offset<2> offset) const
{
  // This function finds the nearest valid pixel along a row or column (specified by 'offset') of an image.
  itk::Index<2> currentPixel = pixel;

  //itk::Size<2> size = this->FullImage->GetLargestPossibleRegion().GetSize();

  // Step forward
  unsigned int pixelCounter = 0;
  while(this->FullImage->GetLargestPossibleRegion().IsInside(currentPixel+offset))
    {
    currentPixel += offset;
    //std::cout << "currentPixel: " << currentPixel << " Valid? "
    //          << this->FullImage->GetPixel(currentPixel).Valid << std::endl;
    pixelCounter++;
    if(this->FullImage->GetPixel(currentPixel).IsValid())
      {
      //std::cout << "Closest valid pixel to " << pixel << " is " << currentPixel << std::endl;
      return currentPixel;
      }
    }

  //std::cout << "Searched " << pixelCounter << " pixels." << std::endl;

  //std::cout << "Offset was " << offset << std::endl;
  offset[0] *= -1;
  offset[1] *= -1;
  //std::cout << "Offset is now " << offset << std::endl;

  currentPixel = pixel; // reset to function input

  // Step backward
  // This is EXACTLY the same loop, but the offset has switched directions
  while(this->FullImage->GetLargestPossibleRegion().IsInside(currentPixel+offset))
    {
    currentPixel += offset;
    //std::cout << "currentPixel: " << currentPixel << " Valid? "
    //          << this->FullImage->GetPixel(currentPixel).Valid << std::endl;
    pixelCounter++;
    if(this->FullImage->GetPixel(currentPixel).IsValid())
      {
      //std::cout << "Closest valid pixel to " << pixel << " is " << currentPixel << std::endl;
      return currentPixel;
      }
    }

  //std::cout << "Searched " << pixelCounter << " pixels." << std::endl;

  // If we get here, an entire row or column did not have a valid pixel. We cannot deal with images this poor.
  std::stringstream ss;
  ss << "No neighbor was found for " << pixel << " in the " << offset << " direction. This should never happen!";
  throw std::runtime_error(ss.str());
  return currentPixel; // So the compiler doesn't complain that all paths don't return a value
}

PTXImage PTXImage::Downsample(const unsigned int factor) const
{
  typedef itk::ShrinkImageFilter <FullImageType, FullImageType> ShrinkImageFilterType;
  ShrinkImageFilterType::Pointer shrinkFilter = ShrinkImageFilterType::New();
  shrinkFilter->SetInput(this->FullImage);
  shrinkFilter->SetShrinkFactor(0, factor);
  shrinkFilter->SetShrinkFactor(1, factor);
  shrinkFilter->Update();

  PTXImage output;
  ITKHelpers::DeepCopy(shrinkFilter->GetOutput(), output.FullImage.GetPointer());

  return output;
}

void PTXImage::WritePTX(const FilePrefix& filePrefix) const
{
  std::stringstream ss;
  ss << filePrefix.prefix << ".ptx";
  std::ofstream fout(ss.str().c_str());

  unsigned int numberOfThetaPoints = this->FullImage->GetLargestPossibleRegion().GetSize()[0];
  unsigned int numberOfPhiPoints = this->FullImage->GetLargestPossibleRegion().GetSize()[1];

  std::cout << "Writing " << ss.str().c_str() << " PTX with: " << std::endl;
  std::cout << "Theta points: " << numberOfThetaPoints << std::endl;
  std::cout << "Phi points: " << numberOfPhiPoints << std::endl;

  fout << numberOfThetaPoints << std::endl;
  fout << numberOfPhiPoints << std::endl;
  fout << "0 0 0" << std::endl
       << "1 0 0" << std::endl
       << "0 1 0" << std::endl
       << "0 0 1" << std::endl
       << "1 0 0 0" << std::endl
       << "0 1 0 0" << std::endl
       << "0 0 1 0" << std::endl
       << "0 0 0 1" << std::endl;

  unsigned int numberOfInvalidPixels = 0;
  unsigned int numberOfZeroPixels = 0;

  for(unsigned int theta = 0; theta < numberOfThetaPoints; theta++)
    {
    for(unsigned int phi = 0; phi < numberOfPhiPoints; phi++)
      {
      itk::Index<2> pixelIndex;
      pixelIndex[0] = theta;
      pixelIndex[1] = phi;

      PTXPixel pixel = this->FullImage->GetPixel(pixelIndex);
      if(pixel.Valid) // Note: this must not be IsValid, as we are testing the IsZero condition separately
        {
        // If the pixel was originally invalid, but has been marked as valid,
        // write it with intensity 1.0 instead of 0.5 so it will be read in
        // as valid in the future.
        if(pixel.IsZero())
          {
          numberOfZeroPixels++;
          fout << pixel.X << " " << pixel.Y << " " << pixel.Z << " " << 1.0
            << " " << static_cast<int>(pixel.R) << " " << static_cast<int>(pixel.G)
            << " " << static_cast<int>(pixel.B) << std::endl;
          }
        else // Write the point as-is
          {
          fout << pixel.X << " " << pixel.Y << " " << pixel.Z << " " << pixel.Intensity
            << " " << static_cast<int>(pixel.R) << " " << static_cast<int>(pixel.G)
            << " " << static_cast<int>(pixel.B) << std::endl;
          }
        }
      else
        {
        numberOfInvalidPixels++;
        /*
        std::cout << "Invalid pixel would have been written as: "
                  << pixel.X << " " << pixel.Y << " " << pixel.Z << " " << pixel.Intensity
                  << " " << static_cast<int>(pixel.R) << " " << static_cast<int>(pixel.G)
                  << " " << static_cast<int>(pixel.B) << std::endl;
        */
        fout << "0 0 0 0.5 0 0 0" << std::endl;
        }
      }
    }

  std::cout << "WritePTX(): Wrote " << numberOfInvalidPixels << " invalid pixels." << std::endl;
  // Note: this will not usually match the CountNumberOfZeroPixels() because it is only pixels that were zero AND valid.
  std::cout << "WritePTX(): Wrote " << numberOfZeroPixels << " valid but zero pixels." << std::endl;

  fout.close();
}

void PTXImage::Crop(const itk::ImageRegion<2>& region)
{
  typedef itk::RegionOfInterestImageFilter< FullImageType, FullImageType > RegionOfInterestImageFilterType;
  RegionOfInterestImageFilterType::Pointer regionOfInterestImageFilter = RegionOfInterestImageFilterType::New();
  regionOfInterestImageFilter->SetRegionOfInterest(region);
  regionOfInterestImageFilter->SetInput(this->FullImage);
  regionOfInterestImageFilter->Update();

  ITKHelpers::DeepCopy(regionOfInterestImageFilter->GetOutput(), this->FullImage.GetPointer());
}

void PTXImage::ComputeWeightedDepthLaplacian(const std::string& filename) const
{
  FloatImageType::Pointer output = FloatImageType::New();
  output->SetRegions(this->FullImage->GetLargestPossibleRegion());
  output->Allocate();

  itk::Size<2> radius;
  radius.Fill(1);
  itk::ConstNeighborhoodIterator<FullImageType> imageIterator(radius, this->FullImage, this->FullImage->GetLargestPossibleRegion());

  FloatImageType::Pointer depthImage = FloatImageType::New();
  CreateDepthImage(depthImage);

  itk::ImageRegionIterator<FloatImageType> outputIterator(output, output->GetLargestPossibleRegion());

  while(!imageIterator.IsAtEnd())
    {
    if(!imageIterator.GetCenterPixel().IsValid())
      {
      ++imageIterator;
      ++outputIterator;
      continue;
      }

    std::vector<itk::Offset<2> > offsets;

    itk::Offset<2> left;
    left[0] = -1;
    left[1] = 0;
    offsets.push_back(left);

    itk::Offset<2> right;
    right[0] = 1;
    right[1] = 0;
    offsets.push_back(right);

    itk::Offset<2> top;
    top[0] = 0;
    top[1] = -1;
    offsets.push_back(top);

    itk::Offset<2> bottom;
    bottom[0] = 0;
    bottom[1] = 1;
    offsets.push_back(bottom);

    std::vector<float> inverseDistances;
    std::vector<float> depths;
    for(unsigned int i = 0; i < offsets.size(); i++)
      {
      if(imageIterator.GetPixel(left).IsValid())
        {
        float dist = DistanceBetweenPoints(imageIterator.GetCenterPixel(), imageIterator.GetPixel(offsets[i]));
        inverseDistances.push_back(1.0/dist);
        depths.push_back(depthImage->GetPixel(imageIterator.GetIndex() + offsets[i]));
        }
      }

    // Normalize inverseDistances
    float sum = 0.0;
    for(unsigned int i = 0; i < inverseDistances.size(); i++)
      {
      sum += inverseDistances[i];
      }
    for(unsigned int i = 0; i < inverseDistances.size(); i++)
      {
      inverseDistances[i] /= sum;
      }
    for(unsigned int i = 0; i < inverseDistances.size(); i++)
      {
      inverseDistances[i] *= static_cast<float>(inverseDistances[i]); // The number of pixels used in the Laplacian computation
      }

    // Compute weighted Laplacian
    float total = 0.0;
    total += depths.size(); // The value of the center pixel in the kernel is the number of neighbors used in the comparison
    for(unsigned int i = 0; i < depths.size(); i++)
      {
      total += inverseDistances[i] * depths[i];
      }

    ++imageIterator;
    ++outputIterator;
    }

  typedef  itk::ImageFileWriter<FloatImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(filename);
  writer->SetInput(output);
  writer->Update();
}

PTXImage::VectorType PTXImage::GetPrincipalAxis() const
{
  // Get center pixel
  itk::Index<2> centerIndex;
  centerIndex[0] = this->FullImage->GetLargestPossibleRegion().GetSize()[0] / 2;
  centerIndex[1] = this->FullImage->GetLargestPossibleRegion().GetSize()[1] / 2;

  FullImageType::PixelType centerPixel = this->FullImage->GetPixel(centerIndex);

  // Assuming the scanner is at the origin, the direction is simply the coordinate
  typedef itk::CovariantVector<double, 3> VectorType;
  VectorType v;
  v[0] = centerPixel.X;
  v[1] = centerPixel.Y;
  v[2] = centerPixel.Z;

  v.Normalize();

  return v;
}

float PTXImage::DistanceBetweenPoints(const PTXPixel& a, const PTXPixel& b) const
{
  if(!a.IsValid() || !b.IsValid())
    {
    std::cerr << "Warning: Comparing distance between pixels which are not both valid." << std::endl;
    return 0.0;
    }

  itk::Point<float, 3> p0;
  p0[0] = a.X;
  p0[1] = a.Y;
  p0[2] = a.Z;

  itk::Point<float, 3> p1;
  p1[0] = b.X;
  p1[1] = b.Y;
  p1[2] = b.Z;

  return p0.EuclideanDistanceTo(p1);
}

void PTXImage::WriteProjectionPlane(const std::string& filename) const
{
  // As always, we assume the scanner is at the origin

  VectorType principalAxis = GetPrincipalAxis();

  double origin[3] = {0,0,0};
  vtkSmartPointer<vtkPlaneSource> planeSource =
    vtkSmartPointer<vtkPlaneSource>::New();
  planeSource->SetCenter(origin);
  planeSource->SetNormal(principalAxis[0], principalAxis[1], principalAxis[2]);
  //planeSource->SetNormal(principalAxis.GetDataPointer()); // only works if the vector is double
  planeSource->Update();

  vtkSmartPointer<vtkLineSource> lineSource =
    vtkSmartPointer<vtkLineSource>::New();
  lineSource->SetPoint1(origin);
  lineSource->SetPoint2(principalAxis[0], principalAxis[1], principalAxis[2]);
  lineSource->Update();

  vtkSmartPointer<vtkAppendPolyData> appendFilter =
    vtkSmartPointer<vtkAppendPolyData>::New();
  appendFilter->AddInputConnection(planeSource->GetOutputPort());
  appendFilter->AddInputConnection(lineSource->GetOutputPort());
  appendFilter->Update();

  vtkSmartPointer<vtkXMLPolyDataWriter> writer =
    vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName(filename.c_str());
  writer->SetInputConnection(appendFilter->GetOutputPort());
  writer->Write();
}

void PTXImage::SetSize(const itk::ImageRegion<2>& region)
{
  this->FullImage->SetRegions(region);
  this->FullImage->Allocate();
  this->FullImage->FillBuffer(PTXPixel());
}

PTXImage PTXImage::OrthogonalProjection(const VectorType& axis) const
{
  //WriteProjectionPlane("projectionPlane.vtp"); // for debugging

  // Create a plane through the origin "aimed" at the scanned points
  vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();
  plane->SetOrigin(0.0, 0.0, 0.0);
  plane->SetNormal(axis[0], axis[1], axis[2]);

  itk::ImageRegionConstIterator<FullImageType> imageIterator(this->FullImage, this->FullImage->GetLargestPossibleRegion());

  vtkSmartPointer<vtkPoints> projectedPoints = vtkSmartPointer<vtkPoints>::New();

  /*
  while(!imageIterator.IsAtEnd())
    {
    FullImageType::PixelType pixel = imageIterator.Get();

    double point[3] = {pixel.X, pixel.Y, pixel.Z};
    double projected[3];
    plane->ProjectPoint(point, projected);

    projectedPoints->InsertNextPoint(projected);

    ++imageIterator;
    }
  */

  // Project the points onto the plane

  vtkSmartPointer<vtkPolyData> pointCloud = vtkSmartPointer<vtkPolyData>::New();
  CreatePointCloud(pointCloud);

  for(unsigned int i = 0; i < static_cast<unsigned int>(pointCloud->GetNumberOfPoints()); i++)
    {
    double point[3];
    pointCloud->GetPoint(i, point);
    double projected[3];
    plane->ProjectPoint(point, projected);

    projectedPoints->InsertNextPoint(projected);

    }

  // Create a polydata of the projected points
  vtkSmartPointer<vtkPolyData> projectedPointsPolydata = vtkSmartPointer<vtkPolyData>::New();
  projectedPointsPolydata->ShallowCopy(pointCloud); // copy the colors and depths
  projectedPointsPolydata->SetPoints(projectedPoints); // replace the point coordinates
  //projectedPointsPolydata->GetPointData()->SetScalars(pointCloud->GetPointData()->GetScalars());

  VTKHelpers::WritePolyData(projectedPointsPolydata, "projectedPoints.vtp");

  //// Determine a coordinate system for the new image ////

  // Compute the projection of the vector pointing to the top-center pixel on the projection plane

  // Get top center pixel
  itk::Index<2> topCenterIndex = FindValidTopCenterPixel();
  std::cout << "topCenterIndex: " << topCenterIndex << std::endl;
  std::cout << "image dims: " << this->FullImage->GetLargestPossibleRegion().GetSize() << std::endl;

  FullImageType::PixelType topCenterPixel = this->FullImage->GetPixel(topCenterIndex);
  std::cout << "topCenterPixel depth: " << topCenterPixel.GetDepth() << std::endl;;
  std::cout << "topCenterPixel: " << topCenterPixel << std::endl;;

  {
  // for debugging only
  topCenterPixel.R = 255;
  topCenterPixel.G = 0;
  topCenterPixel.B = 0;
  this->FullImage->SetPixel(topCenterIndex, topCenterPixel);

  FilePrefix prefix("topCenterColoredRed");
  WriteRGBImage(prefix);
  }

  // Assuming the scanner is at the origin, the direction is simply the normalized coordinate (x-0, y-0, z-0).
  typedef itk::CovariantVector<double, 3> VectorType;
  VectorType topCenterVector;
  topCenterVector[0] = topCenterPixel.X;
  topCenterVector[1] = topCenterPixel.Y;
  topCenterVector[2] = topCenterPixel.Z;
  topCenterVector.Normalize();

  double v[3];
  v[0] = topCenterPixel.X;
  v[1] = topCenterPixel.Y;
  v[2] = topCenterPixel.Z;
  std::cout << "v: " << v[0] << " " << v[1] << " " << v[2] << std::endl;

  double up[3];
  plane->ProjectVector(v, up);
  std::cout << "up: " << up[0] << " " << up[1] << " " << up[2] << std::endl;

  // We will eventually create an image with orientation described by 'up' and the cross product of 'up' and 'principalAxis'

  double horizontalAxis[3];
  double principalAxisArray[3] = {axis[0], axis[1], axis[2]};
  vtkMath::Cross(principalAxisArray, up, horizontalAxis);

  // Create the standard frame
  float standardFrameOrigin[3] = {0,0,0};
  float standardFrameXDirection[3] = {1,0,0};
  float standardFrameYDirection[3] = {0,1,0};
  float standardFrameZDirection[3] = {0,0,1};
  Frame standardFrame(standardFrameOrigin, standardFrameXDirection, standardFrameYDirection, standardFrameZDirection);
  standardFrame.Write("standardFrame.vtp");

  // Get center pixel
  itk::Index<2> centerIndex = FindValidCenterPixel();

  //centerIndex[0] = this->FullImage->GetLargestPossibleRegion().GetSize()[0] / 2;
  //centerIndex[1] = this->FullImage->GetLargestPossibleRegion().GetSize()[1] / 2;

  FullImageType::PixelType centerPixel = this->FullImage->GetPixel(centerIndex);

  {
  // for debugging only
  centerPixel.R = 0;
  centerPixel.G = 255;
  centerPixel.B = 0;
  this->FullImage->SetPixel(centerIndex, centerPixel);
  FilePrefix prefix("centerColoredGreen");
  WriteRGBImage(prefix);
  }

  if(!centerPixel.IsValid())
    {
    throw std::runtime_error("CenterPixel is not valid!");
    
    }
  double centerPixelArray[3] = {centerPixel.X, centerPixel.Y, centerPixel.Z};

  double projectedCenter[3];
  plane->ProjectPoint(centerPixelArray, projectedCenter);

  //float originalFrameOrigin[3] = {projectedCenter[0], projectedCenter[1], projectedCenter[2]};
  float originalFrameOrigin[3] = {0,0,0};
  //float originalFrameXDirection[3] = {horizontalAxis[0], horizontalAxis[1], horizontalAxis[2]};
  float originalFrameXDirection[3] = {static_cast<float>(-horizontalAxis[0]), static_cast<float>(-horizontalAxis[1]), static_cast<float>(-horizontalAxis[2])}; // casts to float are needed to avoid "error: narrowing conversion of double to float inside {}"
  float originalFrameYDirection[3] = {static_cast<float>(up[0]), static_cast<float>(up[1]), static_cast<float>(up[2])};
  float originalFrameZDirection[3] = {static_cast<float>(principalAxisArray[0]), static_cast<float>(principalAxisArray[1]), static_cast<float>(principalAxisArray[2])};
  Frame originalFrame(originalFrameOrigin, originalFrameXDirection, originalFrameYDirection, originalFrameZDirection);
  originalFrame.Write("originalFrame.vtp");

  vtkSmartPointer<vtkTransform> transform =
    vtkSmartPointer<vtkTransform>::New();
  originalFrame.AlignTo(standardFrame, transform);

  vtkSmartPointer<vtkTransformFilter> transformFilter =
    vtkSmartPointer<vtkTransformFilter>::New();
  transformFilter->SetInputData(projectedPointsPolydata);
  transformFilter->SetTransform(transform);
  transformFilter->Update();

  VTKHelpers::WritePolyData(vtkPolyData::SafeDownCast(transformFilter->GetOutput()), "transformedProjectedPoints.vtp");

  // Compute the corner of the image
  double bounds[6];
  transformFilter->GetOutput()->GetBounds(bounds);

  vtkSmartPointer<vtkTransform> translation =
    vtkSmartPointer<vtkTransform>::New();
  float xmin = bounds[0];
  float ymin = bounds[2];
  translation->Translate(-xmin, -ymin, 0);

  vtkSmartPointer<vtkTransformFilter> shiftTransformFilter =
    vtkSmartPointer<vtkTransformFilter>::New();
  shiftTransformFilter->SetInputConnection(transformFilter->GetOutputPort());
  shiftTransformFilter->SetTransform(translation);
  shiftTransformFilter->Update();

  VTKHelpers::WritePolyData(vtkPolyData::SafeDownCast(shiftTransformFilter->GetOutput()), "shiftedTransformedProjectedPoints.vtp");

  // Create an image into which to project the points
  PTXImage orthoPTX;
  orthoPTX.SetSize(this->FullImage->GetLargestPossibleRegion());

  double newBounds[6];
  shiftTransformFilter->GetOutput()->GetBounds(newBounds);

  float binSize[2];
  binSize[0] = newBounds[1] / static_cast<float>(orthoPTX.GetSize()[0]); // xmax / width
  binSize[1] = newBounds[3] / static_cast<float>(orthoPTX.GetSize()[1]); // ymax / height

  // Create an image to track the depths of previously projected pixels
  DepthImageType::Pointer depthImage = DepthImageType::New();
  depthImage->SetRegions(this->FullImage->GetLargestPossibleRegion());
  depthImage->Allocate();
  depthImage->FillBuffer(itk::NumericTraits<float>::max());

  for(unsigned int i = 0; i < static_cast<unsigned int>(shiftTransformFilter->GetOutput()->GetNumberOfPoints()); i++)
    {
    // Divide the 3D coordinates by the number of pixels (aka bins) to get the image coordinate
    double p[3];
    shiftTransformFilter->GetOutput()->GetPoint(i,p);
    itk::Index<2> index;
    index[0] = floor(p[0] / binSize[0]);
    index[1] = floor(p[1] / binSize[1]);
    //std::cout << "Index: " << index << std::endl;

    unsigned char pointColor[3];
    vtkUnsignedCharArray::SafeDownCast(shiftTransformFilter->GetOutput()->GetPointData()->GetScalars())->GetTupleValue(i, pointColor);

    // Get the index of the pixel in the original image from which this point came from
    int originalPixel[2];
    vtkIntArray::SafeDownCast(shiftTransformFilter->GetOutput()->
            GetPointData()->GetArray("OriginalPixel"))->GetTupleValue(i, originalPixel);

    itk::Index<2> originalPixelIndex;
    originalPixelIndex[0] = originalPixel[0];
    originalPixelIndex[1] = originalPixel[1];

    // Create the new PTX point
    PTXPixel newPixel;
    newPixel.Valid = true;
    newPixel.R = pointColor[0];
    newPixel.G = pointColor[1];
    newPixel.B = pointColor[2];
    //newPixel.X = p[0];
    //newPixel.Y = p[1];
    //newPixel.Z = p[2];
    newPixel.X = FullImage->GetPixel(originalPixelIndex).X;
    newPixel.Y = FullImage->GetPixel(originalPixelIndex).Y;
    newPixel.Z = FullImage->GetPixel(originalPixelIndex).Z;

    // Ensure the projected point is the closest to the camera to be projected into this bin
    if(depthImage->GetPixel(index) < FullImage->GetPixel(originalPixelIndex).GetDepth())
      {
      // do nothing, there is already a closer pixel projected
#ifdef DEBUGMODE
      std::cout << "Trying to project pixel from depth " << FullImage->GetPixel(originalPixelIndex).GetDepth()
                << " but there is already a projected pixel from depth " << depthImage->GetPixel(index) << std::endl;
#endif
      }
    else
      {
      // Project the point into this bin
      orthoPTX.SetPixel(index, newPixel);

      // Update the "depth buffer" to track previously projected pixels
      depthImage->SetPixel(index, FullImage->GetPixel(originalPixelIndex).GetDepth());
      }
    }

  return orthoPTX;
}

void PTXImage::SetPixel(const itk::Index<2>& index, const PTXPixel& pixel)
{
  this->FullImage->SetPixel(index, pixel);
}

PTXImage::FullImageType::Pointer PTXImage::GetFullImage() const
{
  return this->FullImage;
}

itk::Size<2> PTXImage::GetSize() const
{
  return this->FullImage->GetLargestPossibleRegion().GetSize();
}

itk::Index<2> PTXImage::FindValidTopCenterPixel() const
{
  itk::Index<2> topCenterIndex;
  topCenterIndex[0] = this->FullImage->GetLargestPossibleRegion().GetSize()[0] / 2;
  topCenterIndex[1] = this->FullImage->GetLargestPossibleRegion().GetSize()[1] - 1;

  if(this->FullImage->GetPixel(topCenterIndex).IsValid())
    {
    return topCenterIndex;
    }

  topCenterIndex[1] -= 1;
  if(this->FullImage->GetPixel(topCenterIndex).IsValid())
    {
    return topCenterIndex;
    }

  topCenterIndex[1] += 1;
  topCenterIndex[0] += 1;
  if(this->FullImage->GetPixel(topCenterIndex).IsValid())
    {
    return topCenterIndex;
    }

  topCenterIndex[0] -= 2;
  if(this->FullImage->GetPixel(topCenterIndex).IsValid())
    {
    return topCenterIndex;
    }

  // if we get to here, all pixels around the top center pixel are invalid
  std::cerr << "All pixels around the top center pixel are invalid!" << std::endl;
  exit(-1);
  return topCenterIndex;
}


itk::Index<2> PTXImage::FindValidCenterPixel() const
{
  itk::Index<2> centerIndex;
  centerIndex[0] = this->FullImage->GetLargestPossibleRegion().GetSize()[0] / 2;
  centerIndex[1] = this->FullImage->GetLargestPossibleRegion().GetSize()[1] / 2;

  if(this->FullImage->GetPixel(centerIndex).IsValid())
    {
    return centerIndex;
    }

  itk::Offset<2> offset;
  offset[0] = 1;
  offset[1] = 0;

  if(this->FullImage->GetPixel(centerIndex + offset).IsValid())
    {
    return centerIndex + offset;
    }

  offset[0] = -1;
  offset[1] = 0;

  if(this->FullImage->GetPixel(centerIndex + offset).IsValid())
    {
    return centerIndex + offset;
    }

  offset[0] = 0;
  offset[1] = 1;
  if(this->FullImage->GetPixel(centerIndex + offset).IsValid())
    {
    return centerIndex + offset;
    }

  offset[0] = 0;
  offset[1] = -1;
  if(this->FullImage->GetPixel(centerIndex + offset).IsValid())
    {
    return centerIndex + offset;
    }

  // if we get to here, all pixels around the center pixel are invalid
  std::cerr << "All pixels around the center pixel are invalid!" << std::endl;
  exit(-1);
  return centerIndex;
}

unsigned int PTXImage::GetHeight() const
{
  return this->GetSize()[1];
}

unsigned int PTXImage::GetWidth() const
{
  return this->GetSize()[0];
}

void PTXImage::Backup()
{
  ITKHelpers::DeepCopy(FullImage.GetPointer(), OriginalFullImage.GetPointer());
}

itk::ImageRegion<2> PTXImage::GetFullRegion() const
{
  return this->FullImage->GetLargestPossibleRegion();
}

// void PTXImage::GetMesh(vtkPolyData* const output)
// {
//   output->DeepCopy(this->Mesh);
// }

vtkPolyData* PTXImage::GetMesh() const
{
  return this->Mesh;
}

void PTXImage::ComputeMesh(const float maxMeshEdgeLength)
{
  // Create a grid of theta/phi coordinates (keeps only connectivity, not geometry)
  vtkSmartPointer<vtkPoints> points2D = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkPoints> points3D = vtkSmartPointer<vtkPoints>::New();

  itk::ImageRegionConstIterator<FullImageType> fullImageIterator(this->FullImage, this->FullImage->GetLargestPossibleRegion());

  while(!fullImageIterator.IsAtEnd())
    {
    if(fullImageIterator.Get().IsValid())
      {
      points2D->InsertNextPoint(fullImageIterator.GetIndex()[0], fullImageIterator.GetIndex()[1], 0.0f);
      points3D->InsertNextPoint(fullImageIterator.Get().X, fullImageIterator.Get().Y, fullImageIterator.Get().Z);
      }
    ++fullImageIterator;
    }

  std::cout << "There are " << points2D->GetNumberOfPoints() << " 2D points and " << points3D->GetNumberOfPoints() << " 3D points." << std::endl;
  
  // Add the 2d grid points to a polydata object
  vtkSmartPointer<vtkPolyData> polydata2d = vtkSmartPointer<vtkPolyData>::New();
  polydata2d->SetPoints(points2D);

  // Triangulate the grid points
  vtkSmartPointer<vtkDelaunay2D> delaunay = vtkSmartPointer<vtkDelaunay2D>::New();
  delaunay->SetInputData(polydata2d);
  delaunay->Update();

  std::cout << "Finished delaunay." << std::endl;
  
  // Get the resulting triangles from the triangulation
  vtkCellArray* cells = delaunay->GetOutput()->GetPolys();

  // Create the 3d triangle array
  vtkSmartPointer<vtkCellArray> Triangles3D = vtkSmartPointer<vtkCellArray>::New();

  // Initialize some variables
  vtkIdType npts; // the number of points in a cell
  vtkIdType* pts; //indexes to the points

  // Go through all the triangles of the Delaunay triangulation and add them to the 3d polydata if they are shorter than MaxMeshEdgeLength
  cells->InitTraversal();
  while (cells->GetNextCell(npts,pts))
    {
    // Get the 3 points of the current triangle
    double p0[3];
    double p1[3];
    double p2[3];

    points3D->GetPoint(pts[0], p0);
    points3D->GetPoint(pts[1], p1);
    points3D->GetPoint(pts[2], p2);

    // Throw away triangles that are bigger than a threshold
    double edgeLength = vtkMath::Distance2BetweenPoints(p0, p1);
    if(edgeLength > maxMeshEdgeLength)
      {
      continue;
      }

    edgeLength = vtkMath::Distance2BetweenPoints(p1, p2);
    if(edgeLength > maxMeshEdgeLength)
      {
      continue;
      }

    edgeLength = vtkMath::Distance2BetweenPoints(p0, p2);
    if(edgeLength > maxMeshEdgeLength)
      {
      continue;
      }

    // Add the triangle to the 3d polydata
    vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();

    triangle->GetPointIds()->SetId(0,pts[0]);
    triangle->GetPointIds()->SetId(1,pts[1]);
    triangle->GetPointIds()->SetId(2,pts[2]);
    Triangles3D->InsertNextCell(triangle);

    }//end while

  std::cout << "Finished adding triangles." << std::endl;
  
  // Save the 3d triangles in the output polydata
  vtkSmartPointer<vtkPolyData> polydata3d = vtkSmartPointer<vtkPolyData>::New();
  polydata3d->SetPoints(points3D);
  this->Mesh->DeepCopy(polydata3d);
  this->Mesh->SetPolys(Triangles3D);

}

void PTXImage::OutputInfo()
{
  std::cout << "Invalid points: " << CountInvalidPoints() << std::endl;
  std::cout << "Valid points: " << CountValidPoints() << std::endl;
  std::cout << "Zero points: " << CountZeroPoints() << std::endl;
}

void PTXImage::PrintCoordinate(const itk::Index<2>& index) const
{
  std::cout << "X: " << this->FullImage->GetPixel(index).X << " y: " <<
                        this->FullImage->GetPixel(index).Y << " z: " << 
                        this->FullImage->GetPixel(index).Z << std::endl;
}
