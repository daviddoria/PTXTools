#include "PTXImage.h"

// STL
#include <sstream>

// Helper functions
static PTXImage CreatePTX();

// Test functions
static void TestCountInvalidPoints();
static void TestCountValidPoints();
static void TestCountZeroPoints();

static void TestGetPhi();
static void TestGetTheta();

static void TestMinPhi();
static void TestMaxPhi();
static void TestMinTheta();
static void TestMaxTheta();

static void TestComputeAverageDeltaPhi();
static void TestComputeAverageDeltaTheta();

static void TestApproximateTheta();
static void TestApproximatePhi();

int main()
{
  TestCountInvalidPoints();
  TestCountValidPoints();
  TestCountZeroPoints();

  TestGetPhi();
  TestGetTheta();

  TestComputeAverageDeltaPhi();
  TestComputeAverageDeltaTheta();

  TestApproximateTheta();
  TestApproximatePhi();

  TestMinPhi();
  TestMaxPhi();
  TestMinTheta();
  TestMaxTheta();

  return 0;
}

void TestCountInvalidPoints()
{
  std::cout << "TestCountInvalidPoints()" << std::endl;
  PTXImage ptxImage = CreatePTX();
}

void TestCountValidPoints()
{
  std::cout << "TestCountValidPoints()" << std::endl;
  PTXImage ptxImage = CreatePTX();
}

void TestCountZeroPoints()
{
  std::cout << "TestCountZeroPoints()" << std::endl;
  PTXImage ptxImage = CreatePTX();
}

void TestGetPhi()
{
  std::cout << "TestGetPhi()" << std::endl;
  itk::Index<2> index = {{0,0}};

  PTXImage ptxImage = CreatePTX();

  float phi = ptxImage.GetPhi(index);

  std::cout << "phi: " << phi << std::endl;
}

void TestGetTheta()
{
  std::cout << "TestGetTheta()" << std::endl;
  itk::Index<2> index = {{0,0}};

  PTXImage ptxImage = CreatePTX();

  float theta = ptxImage.GetTheta(index);

  std::cout << "theta: " << theta << std::endl;
}

void TestComputeAverageDeltaPhi()
{
  std::cout << "TestComputeAverageDeltaPhi()" << std::endl;

  PTXImage ptxImage = CreatePTX();

  ptxImage.ComputeAverageDeltaPhi();
  float averageDeltaPhi = ptxImage.GetAverageDeltaPhi();

  std::cout << "averageDeltaPhi: " << averageDeltaPhi << std::endl;
}

void TestComputeAverageDeltaTheta()
{
  std::cout << "TestComputeAverageDeltaTheta()" << std::endl;

  PTXImage ptxImage = CreatePTX();

  ptxImage.ComputeAverageDeltaTheta();
  float averageDeltatheta = ptxImage.GetAverageDeltaTheta();

  std::cout << "averageDeltatheta: " << averageDeltatheta << std::endl;
}

void TestApproximateTheta()
{
  std::cout << "TestApproximateTheta()" << std::endl;

  PTXImage ptxImage = CreatePTX();

  itk::Index<2> index = {{0,0}};

  float approximateTheta = ptxImage.ApproximateTheta(index);

  std::cout << "approximateTheta: " << approximateTheta << std::endl;
}

void TestApproximatePhi()
{
  std::cout << "TestApproximatePhi()" << std::endl;

  PTXImage ptxImage = CreatePTX();

  itk::Index<2> index = {{0,0}};

  float approximatePhi = ptxImage.ApproximatePhi(index);

  std::cout << "approximatePhi: " << approximatePhi << std::endl;
}

void TestMinPhi()
{
  std::cout << "TestMinPhi()" << std::endl;

  PTXImage ptxImage = CreatePTX();

  float minPhi = ptxImage.MinPhi();

  std::cout << "minPhi: " << minPhi << std::endl;
}

void TestMaxPhi()
{
  std::cout << "TestMaxPhi()" << std::endl;

  PTXImage ptxImage = CreatePTX();

  float maxPhi = ptxImage.MaxPhi();

  std::cout << "maxPhi: " << maxPhi << std::endl;
}

void TestMinTheta()
{
  std::cout << "TestMinTheta()" << std::endl;

  PTXImage ptxImage = CreatePTX();

  float minTheta = ptxImage.MinTheta();

  std::cout << "minTheta: " << minTheta << std::endl;
}

void TestMaxTheta()
{
  std::cout << "TestMaxTheta()" << std::endl;

  PTXImage ptxImage = CreatePTX();

  float maxTheta = ptxImage.MaxTheta();

  std::cout << "maxTheta: " << maxTheta << std::endl;
}

// Helper functions
PTXImage CreatePTX()
{
  PTXImage ptxImage;

  // Setup the image
  itk::Size<2> size;
  unsigned int numberOfThetaPoints = 5; // width
  unsigned int numberOfPhiPoints = 3; // height

  size[0] = numberOfThetaPoints;
  size[1] = numberOfPhiPoints;

  itk::Index<2> corner;
  corner.Fill(0);

  itk::ImageRegion<2> region(corner, size);

  ptxImage.GetFullImage()->SetRegions(region);
  ptxImage.GetFullImage()->Allocate();

  itk::ImageRegionIterator<PTXImage::FullImageType> fullImageIterator(ptxImage.GetFullImage(), ptxImage.GetFullImage()->GetLargestPossibleRegion());

  while(!fullImageIterator.IsAtEnd())
  {
    PTXPixel ptxPixel;
    ptxPixel.Valid = true;

    /** The color of the point. */
    ptxPixel.R = 255;
    ptxPixel.G = 255;
    ptxPixel.B = 0;

    /** The Intensity of the reflection at the point. */
    ptxPixel.Intensity = .5;

    /** The Coordinate of the point. */
    ptxPixel.X = fullImageIterator.GetIndex()[0];
    ptxPixel.Y = fullImageIterator.GetIndex()[1];
    ptxPixel.Z = 10;

    fullImageIterator.Set(ptxPixel);
    ++fullImageIterator;
  }

  return ptxImage;
}
