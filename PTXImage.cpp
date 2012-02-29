#include "PTXImage.h"
#include "Frame.h"
#include "Helpers.h"

// ITK
#include "itkAzimuthElevationToCartesianTransform.h"
#include "itkCompose3DCovariantVectorImageFilter.h"
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
#include <vtkFloatArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkPlaneSource.h>
#include <vtkLineSource.h>
#include <vtkPlane.h>
#include <vtkMath.h>
#include <vtkStructuredGrid.h>
#include <vtkUnsignedCharArray.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLStructuredGridWriter.h>

// STL
#include <fstream>
#include <sstream>
#include <string>
#include <Mask.h>

PTXImage::PTXImage()
{
  // Create the main image
  this->FullImage = FullImageType::New();

  this->Debug = false;
}

void PTXImage::SetDebug(bool value)
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
    std::cout << "Pixel " << pixel << " is not inside FullImage: " << this->FullImage->GetLargestPossibleRegion() << std::endl;
    exit(-1);
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

  Helpers::DeepCopy<FullImageType>(tileImageFilter->GetOutput(), this->FullImage);
}

unsigned int PTXImage::CountValidPoints() const
{
  unsigned int totalPoints = this->FullImage->GetLargestPossibleRegion().GetSize()[0] *
                                this->FullImage->GetLargestPossibleRegion().GetSize()[1];
  return  totalPoints - CountInvalidPoints();
}

unsigned int PTXImage::CountInvalidPoints() const
{
  itk::ImageRegionConstIterator<FullImageType> fullImageIterator(this->FullImage, this->FullImage->GetLargestPossibleRegion());

  unsigned int numberOfInvalidPoints = 0;
  while(!fullImageIterator.IsAtEnd())
    {
    PTXPixel fullPixel = fullImageIterator.Get();

    if(!fullPixel.Valid)
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
  typedef itk::Compose3DCovariantVectorImageFilter<FloatImageType,
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
  typedef itk::Compose3DCovariantVectorImageFilter<FloatImageType,
                              XYZImageType> ComposeCovariantVectorImageFilterType;

  ComposeCovariantVectorImageFilterType::Pointer composeFilter = ComposeCovariantVectorImageFilterType::New();

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

    if(fullPixel.Valid)
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

    if(fullPixel.Valid)
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

void PTXImage::CreateRGBImage(RGBImageType::Pointer image) const
{
  // Setup the image
  image->SetRegions(this->FullImage->GetLargestPossibleRegion());
  image->Allocate();
  image->FillBuffer(itk::NumericTraits<RGBImageType::PixelType>::Zero);

  // Setup iterators
  itk::ImageRegionIterator<RGBImageType> rgbImageIterator(image, image->GetLargestPossibleRegion());
  rgbImageIterator.GoToBegin();

  itk::ImageRegionIterator<FullImageType> fullImageIterator(this->FullImage, this->FullImage->GetLargestPossibleRegion());
  fullImageIterator.GoToBegin();

  // Traverse the full image, extracting the color information and saving it in the RGB image
  while(!rgbImageIterator.IsAtEnd())
    {
    PTXPixel fullPixel = fullImageIterator.Get();

    //itk::CovariantVector<unsigned char, 3> rgbPixel;
    RGBImageType::PixelType rgbPixel;

    if(fullPixel.Valid)
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

    if(pixel.Valid)
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
  vertexGlyphFilter->AddInput(polydata);
  vertexGlyphFilter->Update();

  // Add the ptx image size to the field data
  vtkSmartPointer<vtkIntArray> imageSizeArray = vtkSmartPointer<vtkIntArray>::New();
  int imageSize[2] = {this->GetSize()[0], this->GetSize()[1]};
  imageSizeArray->SetNumberOfComponents(2);
  imageSizeArray->SetName("ImageSize");
  imageSizeArray->InsertNextTupleValue(imageSize);
  polydata->GetFieldData()->AddArray(imageSizeArray);

  pointCloud->ShallowCopy(vertexGlyphFilter->GetOutput());
}


void PTXImage::CreateStructuredGrid(vtkSmartPointer<vtkStructuredGrid> structuredGrid) const
{
  int dimensions[3] = {this->GetWidth(), this->GetHeight(), 1};
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

    if(pixel.Valid)
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
  writer->SetInputConnection(pointCloud->GetProducerPort());
  writer->Write();
}

void PTXImage::WriteStructuredGrid(const std::string& fileName) const
{
  vtkSmartPointer<vtkStructuredGrid> structuredGrid =
    vtkSmartPointer<vtkStructuredGrid>::New();
  CreateStructuredGrid(structuredGrid);

  vtkSmartPointer<vtkXMLStructuredGridWriter> writer = vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
  writer->SetFileName(fileName.c_str());
  writer->SetInputConnection(structuredGrid->GetProducerPort());
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
    if(pixel.Valid)
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

  CountInvalidPoints();

  // Setup iterators
  itk::ImageRegionIterator<FullImageType> imageIterator(this->FullImage, this->FullImage->GetLargestPossibleRegion());
  itk::ImageRegionConstIterator<FloatImageType> newDepthIterator(depthImage, depthImage->GetLargestPossibleRegion());

  itk::Point<float, 3> origin;
  origin.Fill(0);

  while(!imageIterator.IsAtEnd())
    {
    // Get the old point
    PTXPixel pixel = imageIterator.Get();
    itk::Point<float, 3> oldPoint;
    oldPoint[0] = pixel.X;
    oldPoint[1] = pixel.Y;
    oldPoint[2] = pixel.Z;

    // For testing only, compute the old depth
    double oldDepth = origin.EuclideanDistanceTo(oldPoint);
    if(this->Debug)
      {
      std::cout << "Old depth: " << oldDepth << " New depth: " << newDepthIterator.Get() << std::endl;
      }

    // Get the vector from the origin (scanner location) and the old point
    itk::Vector<float, 3> unitVector = oldPoint - origin;

    // Get a unit vector in the direction of the old point
    unitVector.Normalize();

    // Compute the new point from the vector and the new depth
    itk::Point<float, 3> newPoint = origin + unitVector * newDepthIterator.Get();

    if(this->Debug)
      {
      std::cout << "Old point: " << oldPoint << " New point: " << newPoint << std::endl;
      }

    // Save the new point in the PTXPixel
    pixel.X = newPoint[0];
    pixel.Y = newPoint[1];
    pixel.Z = newPoint[2];

    imageIterator.Set(pixel);

    ++imageIterator;
    ++newDepthIterator;
    }

}

void PTXImage::ReplaceRGB(const RGBVectorImageType* const rgb)
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
  itk::ImageRegionIterator<FullImageType> imageIterator(this->FullImage, this->FullImage->GetLargestPossibleRegion());
  itk::ImageRegionConstIterator<RGBVectorImageType> rgbIterator(rgb, rgb->GetLargestPossibleRegion());

  while(!imageIterator.IsAtEnd())
    {
    // Get the old point
    PTXPixel pixel = imageIterator.Get();
    if(pixel.Valid)
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
  itk::ImageRegionIterator<FullImageType> imageIterator(this->FullImage, this->FullImage->GetLargestPossibleRegion());
  itk::ImageRegionConstIterator<RGBImageType> rgbIterator(rgb, rgb->GetLargestPossibleRegion());

  while(!imageIterator.IsAtEnd())
    {
    // Get the old point
    PTXPixel pixel = imageIterator.Get();
    if(pixel.Valid)
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
  itk::ImageRegionIteratorWithIndex<FullImageType> imageIterator(this->FullImage, this->FullImage->GetLargestPossibleRegion());

  while(!imageIterator.IsAtEnd())
    {
    // Get the old point
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

  itk::ImageRegionIteratorWithIndex<FullImageType> imageIterator(this->FullImage, this->FullImage->GetLargestPossibleRegion());

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
//   itk::ImageRegionIterator<FullImageType> fullImageIterator(this->FullImage, this->FullImage->GetLargestPossibleRegion());
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
  itk::ImageRegionIterator<FullImageType> imageIterator(this->FullImage, this->FullImage->GetLargestPossibleRegion());
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

void PTXImage::ReplaceRGBD(const RGBDImageType* const rgbd)
{
  // This function allows the color and depth to be modified externally and the new map applied to the grid

  CountInvalidPoints();

  ComputeAverageDeltaPhi();
  ComputeAverageDeltaTheta();

  // Setup iterators
  itk::ImageRegionIterator<FullImageType> imageIterator(this->FullImage, this->FullImage->GetLargestPossibleRegion());
  itk::ImageRegionConstIterator<RGBDImageType> rgbdIterator(rgbd, rgbd->GetLargestPossibleRegion());

  itk::Point<float, 3> origin;
  origin.Fill(0);

  unsigned int madeValid = 0;

  while(!imageIterator.IsAtEnd())
    {
    // Get the old point
    PTXPixel pixel = imageIterator.Get();

    // Copy the color from the RGBD image
    pixel.R = rgbdIterator.Get()[0];
    pixel.G = rgbdIterator.Get()[1];
    pixel.B = rgbdIterator.Get()[2];

    itk::Point<float, 3> oldPoint;

    // If the point was invalid, construct a vector in its direction incase the depth is now valid
    bool previousValidity = pixel.Valid;

    if(pixel.Valid == false)
      {
      madeValid++;
      pixel.Valid = true;

      oldPoint = ApproximateOldPoint(imageIterator.GetIndex());
      std::cout << "new oldPoint: " << oldPoint << std::endl;

      // For debugging, make new pixels bright green
//       pixel.R = 0;
//       pixel.G = 255;
//       pixel.B = 0;

      // For debugging, skip invalid pixels
      ++imageIterator;
      ++rgbdIterator;
      continue;
      }
    else
      {
      // For debugging, turn off the pixels which are not newly valid
      //pixel.Valid = false; //!!!

      oldPoint[0] = pixel.X;
      oldPoint[1] = pixel.Y;
      oldPoint[2] = pixel.Z;
      }

    // Get the vector from the origin (scanner location) and the old point
    itk::Vector<float, 3> unitVector = oldPoint - origin;

    // Get a unit vector in the direction of the old point
    unitVector.Normalize();

    // Compute the new point from the vector and the new depth
    float multiplier = 1.0;
    if(!previousValidity)
      {
      multiplier = -1.0;
      }
    float depth = multiplier * rgbdIterator.Get()[3]; // This works, but why is the negative necessary?
    itk::Point<float, 3> newPoint = origin + unitVector * depth;

    // Save the new point in the PTXPixel
    pixel.X = newPoint[0];
    pixel.Y = newPoint[1];
    pixel.Z = newPoint[2];

    imageIterator.Set(pixel);
    //std::cout << pixel << std::endl;

    ++imageIterator;
    ++rgbdIterator;
    }

  std::cout << madeValid << " pixels were made valid." << std::endl;
  CountInvalidPoints();
}

itk::Point<float, 3> PTXImage::ApproximateOldPoint(const itk::Index<2>& pixel) const
{
  // This function creates a point unit distance from the origin in the approximate direction the old point would have been aquired.
  float phi = ApproximatePhi(pixel);
  float theta = ApproximateTheta(pixel);

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

      if(!this->FullImage->GetPixel(index).Valid || // the current pixel is invalid
        !this->FullImage->GetPixel(lastIndex).Valid || // the previous pixel was invalid
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
    total += deltaThetas[i];
    }

  if(deltaThetas.size() == 0)
    {
    std::cerr << "There are no valid theta neighbors!" << std::endl;
    exit(-1);
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

      if(!this->FullImage->GetPixel(index).Valid || // the current pixel is invalid
        !this->FullImage->GetPixel(lastIndex).Valid || // the previous pixel is invalid
        lastIndex[0] != index[0] // if we are not on the same column as the last pixel
        )
        {
        // do nothing
        }
      else
        {
        deltaPhis.push_back(GetPhi(index) - GetPhi(lastIndex));
        }

      lastIndex = index;
      }
    }

  float total = 0.0;
  for(unsigned int i = 0; i < deltaPhis.size(); i++)
    {
    total += deltaPhis[i];
    }

  if(deltaPhis.size() == 0)
    {
    std::cerr << "There are no valid phi neighbors!" << std::endl;
    exit(-1);
    }

  this->AverageDeltaPhi = total/static_cast<float>(deltaPhis.size());

  std::cout << "AverageDeltaPhi: " << this->AverageDeltaPhi << std::endl;
}

float PTXImage::GetPhi(const itk::Index<2>& index) const
{
  if(!this->FullImage->GetPixel(index).Valid)
    {
    std::cerr << "Cannot GetPhi on an invalid pixel!" << std::endl;
    exit(-1);
    }

  typedef itk::Point<float, 3> PointType;
  PointType cartesian;
  cartesian[0] = this->FullImage->GetPixel(index).X;
  cartesian[1] = this->FullImage->GetPixel(index).Y;
  cartesian[2] = this->FullImage->GetPixel(index).Z;

  typedef itk::AzimuthElevationToCartesianTransform< float, 3 >
    AzimuthElevationToCartesian;
  AzimuthElevationToCartesian::Pointer azimuthElevation =
    AzimuthElevationToCartesian::New();

  return azimuthElevation->TransformCartesianToAzEl(cartesian)[1];
}

float PTXImage::GetTheta(const itk::Index<2>& index) const
{
  if(!this->FullImage->GetPixel(index).Valid)
    {
    std::cerr << "Cannot GetPhi on an invalid pixel!" << std::endl;
    exit(-1);
    }

  typedef itk::Point<float, 3> PointType;
  PointType cartesian;
  cartesian[0] = this->FullImage->GetPixel(index).X;
  cartesian[1] = this->FullImage->GetPixel(index).Y;
  cartesian[2] = this->FullImage->GetPixel(index).Z;

  typedef itk::AzimuthElevationToCartesianTransform< float, 3 >
    AzimuthElevationToCartesian;
  AzimuthElevationToCartesian::Pointer azimuthElevation =
    AzimuthElevationToCartesian::New();

  return azimuthElevation->TransformCartesianToAzEl(cartesian)[0];
}

float PTXImage::ApproximateTheta(const itk::Index<2>& pixel) const
{
  itk::Offset<2> offset;
  offset[0] = 1;
  offset[1] = 0;

  itk::Index<2> nearestPixel = FindNearestValidPixel(pixel, offset);

  float theta = GetTheta(nearestPixel);

  return theta + (pixel[0] - nearestPixel[0])*this->AverageDeltaTheta;
}

float PTXImage::ApproximatePhi(const itk::Index<2>& pixel) const
{
  itk::Offset<2> offset;
  offset[0] = 0;
  offset[1] = 1;

  itk::Index<2> nearestPixel = FindNearestValidPixel(pixel, offset);

  float phi = GetPhi(nearestPixel);

  return phi + (pixel[1] - nearestPixel[1])*this->AverageDeltaPhi;
}

itk::Index<2> PTXImage::FindNearestValidPixel(const itk::Index<2>& pixel, const itk::Offset<2>& inputOffset) const
{
  // This function finds the nearest valid pixel along a row or column (specified by 'offset') of an image.
  itk::Index<2> currentPixel = pixel;
  itk::Offset<2> offset = inputOffset; // We will need to modify this internally, so we need to create a new object so the input can be const.

  //itk::Size<2> size = this->FullImage->GetLargestPossibleRegion().GetSize();

  // Step forward
  unsigned int pixelCounter = 0;
  while(this->FullImage->GetLargestPossibleRegion().IsInside(currentPixel+offset))
    {
    currentPixel += offset;
    //std::cout << "currentPixel: " << currentPixel << " Valid? " << this->FullImage->GetPixel(currentPixel).Valid << std::endl;
    pixelCounter++;
    if(this->FullImage->GetPixel(currentPixel).Valid)
      {
      std::cout << "Closest valid pixel to " << pixel << " is " << currentPixel << std::endl;
      return currentPixel;
      }
    }

  std::cout << "Searched " << pixelCounter << " pixels." << std::endl;

  std::cout << "Offset was " << offset << std::endl;
  offset[0] *= -1;
  offset[1] *= -1;
  std::cout << "Offset is now " << offset << std::endl;

  currentPixel = pixel; // reset to function input

  // Step backward
  // This is EXACTLY the same loop, but the offset has switched directions
  while(this->FullImage->GetLargestPossibleRegion().IsInside(currentPixel+offset))
    {
    currentPixel += offset;
    //std::cout << "currentPixel: " << currentPixel << " Valid? " << this->FullImage->GetPixel(currentPixel).Valid << std::endl;
    pixelCounter++;
    if(this->FullImage->GetPixel(currentPixel).Valid)
      {
      std::cout << "Closest valid pixel to " << pixel << " is " << currentPixel << std::endl;
      return currentPixel;
      }
    }

  //std::cout << "Searched " << pixelCounter << " pixels." << std::endl;

  // If we get here, an entire row or column did not have a valid pixel. We cannot deal with images this poor.
  std::cerr << "No neighbor was found for " << pixel << " in the " << offset << " direction. This should never happen!" << std::endl;
  exit(-1);
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
  Helpers::DeepCopy<FullImageType>(shrinkFilter->GetOutput(), output.FullImage);

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

  for(unsigned int theta = 0; theta < numberOfThetaPoints; theta++)
    {
    for(unsigned int phi = 0; phi < numberOfPhiPoints; phi++)
      {
      itk::Index<2> pixelIndex;
      pixelIndex[0] = theta;
      pixelIndex[1] = phi;

      PTXPixel pixel = this->FullImage->GetPixel(pixelIndex);
      if(pixel.Valid)
        {
        fout << pixel.X << " " << pixel.Y << " " << pixel.Z << " " << pixel.Intensity
            << " " << static_cast<int>(pixel.R) << " " << static_cast<int>(pixel.G)
            << " " << static_cast<int>(pixel.B) << std::endl;
        }
      else
        {
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

  fout.close();
}

void PTXImage::Crop(const itk::ImageRegion<2>& region)
{
  typedef itk::RegionOfInterestImageFilter< FullImageType, FullImageType > RegionOfInterestImageFilterType;
  RegionOfInterestImageFilterType::Pointer regionOfInterestImageFilter = RegionOfInterestImageFilterType::New();
  regionOfInterestImageFilter->SetRegionOfInterest(region);
  regionOfInterestImageFilter->SetInput(this->FullImage);
  regionOfInterestImageFilter->Update();

  Helpers::DeepCopy<FullImageType>(regionOfInterestImageFilter->GetOutput(), this->FullImage);
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
    if(!imageIterator.GetCenterPixel().Valid)
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
      if(imageIterator.GetPixel(left).Valid)
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
  if(!a.Valid || !b.Valid)
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

  Helpers::OutputPolyData(projectedPointsPolydata, "projectedPoints.vtp");

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

  if(!centerPixel.Valid)
    {
    std::cerr << "CenterPixel is not valid!" << std::endl;
    exit(-1);
    }
  double centerPixelArray[3] = {centerPixel.X, centerPixel.Y, centerPixel.Z};

  double projectedCenter[3];
  plane->ProjectPoint(centerPixelArray, projectedCenter);

  //float originalFrameOrigin[3] = {projectedCenter[0], projectedCenter[1], projectedCenter[2]};
  float originalFrameOrigin[3] = {0,0,0};
  //float originalFrameXDirection[3] = {horizontalAxis[0], horizontalAxis[1], horizontalAxis[2]};
  float originalFrameXDirection[3] = {-horizontalAxis[0], -horizontalAxis[1], -horizontalAxis[2]};
  float originalFrameYDirection[3] = {up[0], up[1], up[2]};
  float originalFrameZDirection[3] = {principalAxisArray[0], principalAxisArray[1], principalAxisArray[2]};
  Frame originalFrame(originalFrameOrigin, originalFrameXDirection, originalFrameYDirection, originalFrameZDirection);
  originalFrame.Write("originalFrame.vtp");

  vtkSmartPointer<vtkTransform> transform =
    vtkSmartPointer<vtkTransform>::New();
  originalFrame.AlignTo(standardFrame, transform);

  vtkSmartPointer<vtkTransformFilter> transformFilter =
    vtkSmartPointer<vtkTransformFilter>::New();
  transformFilter->SetInputConnection(projectedPointsPolydata->GetProducerPort());
  transformFilter->SetTransform(transform);
  transformFilter->Update();

  Helpers::OutputPolyData(vtkPolyData::SafeDownCast(transformFilter->GetOutput()), "transformedProjectedPoints.vtp");

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

  Helpers::OutputPolyData(vtkPolyData::SafeDownCast(shiftTransformFilter->GetOutput()), "shiftedTransformedProjectedPoints.vtp");

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
    vtkIntArray::SafeDownCast(shiftTransformFilter->GetOutput()->GetPointData()->GetArray("OriginalPixel"))->GetTupleValue(i, originalPixel);

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

  if(this->FullImage->GetPixel(topCenterIndex).Valid)
    {
    return topCenterIndex;
    }

  topCenterIndex[1] -= 1;
  if(this->FullImage->GetPixel(topCenterIndex).Valid)
    {
    return topCenterIndex;
    }

  topCenterIndex[1] += 1;
  topCenterIndex[0] += 1;
  if(this->FullImage->GetPixel(topCenterIndex).Valid)
    {
    return topCenterIndex;
    }

  topCenterIndex[0] -= 2;
  if(this->FullImage->GetPixel(topCenterIndex).Valid)
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

  if(this->FullImage->GetPixel(centerIndex).Valid)
    {
    return centerIndex;
    }

  itk::Offset<2> offset;
  offset[0] = 1;
  offset[1] = 0;

  if(this->FullImage->GetPixel(centerIndex + offset).Valid)
    {
    return centerIndex + offset;
    }

  offset[0] = -1;
  offset[1] = 0;

  if(this->FullImage->GetPixel(centerIndex + offset).Valid)
    {
    return centerIndex + offset;
    }

  offset[0] = 0;
  offset[1] = 1;
  if(this->FullImage->GetPixel(centerIndex + offset).Valid)
    {
    return centerIndex + offset;
    }

  offset[0] = 0;
  offset[1] = -1;
  if(this->FullImage->GetPixel(centerIndex + offset).Valid)
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
