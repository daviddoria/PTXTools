
void PTXImage::WritePointCloud(std::string filename)
{
  // Compute angular spacing
  itk::Index<2> cornerPixelIndex;
  cornerPixelIndex.Fill(0);

  PTXPixel corner = this->FullImage->GetPixel(cornerPixelIndex);
  float theta0 = corner.GetTheta();
  float phi0 = corner.GetPhi();

  std::cout << "theta start: " << corner.GetTheta() << std::endl;
  std::cout << "phi start: " << corner.GetPhi() << std::endl;

  itk::Index<2> index1;
  index1[0] = 1;
  index1[1] = 0;
  PTXPixel pixel1 = this->FullImage->GetPixel(index1);
  float theta1 = pixel1.GetTheta();

  itk::Index<2> index2;
  index2[0] = 0;
  index2[1] = 1;
  PTXPixel pixel2 = this->FullImage->GetPixel(index2);
  float phi1 = pixel2.GetPhi();

  float thetaStep = fabs(theta0 - theta1);
  float phiStep = fabs(phi0 - phi1);

  std::cout << "theta step: " << thetaStep << std::endl;
  std::cout << "phi step: " << phiStep << std::endl;

  vtkSmartPointer<vtkPolyData> pointCloud =
    vtkSmartPointer<vtkPolyData>::New();

  vtkSmartPointer<vtkPoints> points =
    vtkSmartPointer<vtkPoints>::New();

  vtkSmartPointer<vtkUnsignedCharArray> colors =
    vtkSmartPointer<vtkUnsignedCharArray>::New();
  colors->SetNumberOfComponents(3);
  colors->SetName("Colors");

  itk::ImageRegionConstIteratorWithIndex<FullImageType> imageIterator(this->FullImage, this->FullImage->GetLargestPossibleRegion());

  while(!imageIterator.IsAtEnd())
    {
    // Get the value of the current pixel
    PTXPixel pixel = imageIterator.Get();
    unsigned char rgb[3];
    rgb[0] = pixel[0];
    rgb[1] = pixel[1];
    rgb[2] = pixel[2];
    colors->InsertNextTupleValue(rgb);

    float theta = thetaStart + imageIterator.GetIndex()[0] * thetaStep;
    float phi = phiStart + imageIterator.GetIndex()[1] * phiStep;
    float rho = pixel[3];

    float x = rho * sin(theta) * cos(phi);
    float y = rho * sin(theta) * sin(phi);
    float z = rho * sin(phi);

    points->InsertNextPoint(x,y,z);

    ++imageIterator;
    }

  vtkSmartPointer<vtkPolyData> polydata =
    vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints(points);
  polydata->GetPointData()->SetScalars(colors);

  vtkSmartPointer<vtkVertexGlyphFilter> vertexGlyphFilter =
    vtkSmartPointer<vtkVertexGlyphFilter>::New();
  vertexGlyphFilter->AddInput(polydata);
  vertexGlyphFilter->Update();

  vtkSmartPointer<vtkXMLPolyDataWriter> writer =
    vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName(outputFilename.c_str());
  writer->SetInputConnection(pointCloud->GetProducerPort());
  writer->Write();
}