#ifndef FRAME_H
#define FRAME_H

#include <vtkAppendPolyData.h>
#include <vtkCellData.h>
#include <vtkLandmarkTransform.h>
#include <vtkLineSource.h>
#include <vtkMath.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkTransform.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkTransformFilter.h>

struct Frame
{
  Frame(float o[3], float x[3], float y[3], float z[3])
  {
    this->SetOrigin(o);
    this->SetXDirection(x);
    this->SetYDirection(y);
    this->SetZDirection(z);

    std::cout << "Origin: " << this->origin[0] << " " << this->origin[1] << " " << this->origin[2] << std::endl;
    std::cout << "xDirection: "<< this->xDirection[0] << " " << this->xDirection[1] << " " << this->xDirection[2] << std::endl;
    std::cout << "yDirection: "<< this->yDirection[0] << " " << this->yDirection[1] << " " << this->yDirection[2] << std::endl;
    std::cout << "zDirection: "<< this->zDirection[0] << " " << this->zDirection[1] << " " << this->zDirection[2] << std::endl;
  }

  void ApplyTransform(vtkTransform* transform, std::string filename)
  {
    vtkSmartPointer<vtkPolyData> polydata =
      vtkSmartPointer<vtkPolyData>::New();
    CreatePolydata(polydata);

    vtkSmartPointer<vtkTransformFilter> transformFilter =
      vtkSmartPointer<vtkTransformFilter>::New();
    transformFilter->SetInputConnection(polydata->GetProducerPort());
    transformFilter->SetTransform(transform);
    transformFilter->Update();

    vtkSmartPointer<vtkXMLPolyDataWriter> writer =
      vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputConnection(transformFilter->GetOutputPort());
    writer->Write();
  }

  void CreatePolydata(vtkPolyData* polydata)
  {
    /*
    // Create end points of the frame only
    vtkSmartPointer<vtkPoints> points =
      vtkSmartPointer<vtkPoints>::New();

    points->InsertNextPoint(this->origin);
    float x[3];
    vtkMath::Add(this->origin, this->xDirection, x);
    points->InsertNextPoint(x);
    float y[3];
    vtkMath::Add(this->origin, this->yDirection, y);
    points->InsertNextPoint(y);
    float z[3];
    vtkMath::Add(this->origin, this->zDirection, z);
    points->InsertNextPoint(z);

    polydata->SetPoints(points);

    vtkSmartPointer<vtkVertexGlyphFilter> vertexGlyphFilter =
      vtkSmartPointer<vtkVertexGlyphFilter>::New();
    vertexGlyphFilter->AddInput(polydata);
    vertexGlyphFilter->Update();

    polydata->ShallowCopy(vertexGlyphFilter->GetOutput());
    */

    // Create lines of the frame
    float x[3];
    vtkMath::Add(this->origin, this->xDirection, x);
    std::cout << "X: " << x[0] << " " << x[1] << " " << x[2] << std::endl;

    float y[3];
    vtkMath::Add(this->origin, this->yDirection, y);
    std::cout << "Y: " << y[0] << " " << y[1] << " " << y[2] << std::endl;

    float z[3];
    vtkMath::Add(this->origin, this->zDirection, z);
    std::cout << "Z: " << z[0] << " " << z[1] << " " << z[2] << std::endl;

    vtkSmartPointer<vtkLineSource> xAxis =
      vtkSmartPointer<vtkLineSource>::New();
    //xAxis->SetPoint1(origin); //need to add a vtkLineSource::SetPoint1(float*)
    xAxis->SetPoint1(origin[0], origin[1], origin[2]);
    xAxis->SetPoint2(x[0], x[1], x[2]);
    xAxis->Update();

    vtkSmartPointer<vtkLineSource> yAxis =
      vtkSmartPointer<vtkLineSource>::New();
    //yAxis->SetPoint1(origin);
    yAxis->SetPoint1(origin[0], origin[1], origin[2]);
    yAxis->SetPoint2(y[0], y[1], y[2]);
    yAxis->Update();

    vtkSmartPointer<vtkLineSource> zAxis =
      vtkSmartPointer<vtkLineSource>::New();
    zAxis->SetPoint1(origin[0], origin[1], origin[2]);
    zAxis->SetPoint2(z[0], z[1], z[2]);
    zAxis->Update();

    vtkSmartPointer<vtkAppendPolyData> appendFilter =
      vtkSmartPointer<vtkAppendPolyData>::New();
    appendFilter->AddInputConnection(xAxis->GetOutputPort());
    appendFilter->AddInputConnection(yAxis->GetOutputPort());
    appendFilter->AddInputConnection(zAxis->GetOutputPort());
    appendFilter->Update();

    polydata->ShallowCopy(appendFilter->GetOutput());

    // Setup two colors - one for each line
    unsigned char red[3] = {255, 0, 0};
    unsigned char green[3] = {0, 255, 0};
    unsigned char blue[3] = {0, 0, 255};

    // Setup the colors array
    vtkSmartPointer<vtkUnsignedCharArray> colors =
      vtkSmartPointer<vtkUnsignedCharArray>::New();
    colors->SetNumberOfComponents(3);
    colors->SetName("Colors");

    // Add the colors we created to the colors array
    colors->InsertNextTupleValue(red);
    colors->InsertNextTupleValue(green);
    colors->InsertNextTupleValue(blue);

    polydata->GetCellData()->SetScalars(colors);
  }

  void Write(std::string filename)
  {
    vtkSmartPointer<vtkPolyData> polydata =
      vtkSmartPointer<vtkPolyData>::New();
    CreatePolydata(polydata);

    vtkSmartPointer<vtkXMLPolyDataWriter> writer =
      vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputConnection(polydata->GetProducerPort());
    writer->Write();
  }

  float origin[3];
  float xDirection[3];
  float yDirection[3];
  float zDirection[3];

  void SetOrigin(float o[3])
  {
    this->origin[0] = o[0];
    this->origin[1] = o[1];
    this->origin[2] = o[2];
  }

  void SetXDirection(float direction[3])
  {
    vtkMath::Normalize(direction);
    this->xDirection[0] = direction[0];
    this->xDirection[1] = direction[1];
    this->xDirection[2] = direction[2];
  }

  void SetYDirection(float direction[3])
  {
    vtkMath::Normalize(direction);
    this->yDirection[0] = direction[0];
    this->yDirection[1] = direction[1];
    this->yDirection[2] = direction[2];
  }

  void SetZDirection(float direction[3])
  {
    vtkMath::Normalize(direction);
    this->zDirection[0] = direction[0];
    this->zDirection[1] = direction[1];
    this->zDirection[2] = direction[2];
  }

  void AlignTo(const Frame targetFrame, vtkTransform* transform)
  {
    // This function aligns the 'this' frame with targetFrame and returns the transform that was used in 'transform'

    vtkSmartPointer<vtkLandmarkTransform> landmarkTransform =
      vtkSmartPointer<vtkLandmarkTransform>::New();

    // Setup source points
    vtkSmartPointer<vtkPoints> sourcePoints =
      vtkSmartPointer<vtkPoints>::New();

    sourcePoints->InsertNextPoint(this->origin);
    float sourceX[3];
    vtkMath::Add(this->origin, this->xDirection, sourceX);
    sourcePoints->InsertNextPoint(sourceX);
    float sourceY[3];
    vtkMath::Add(this->origin, this->yDirection, sourceY);
    sourcePoints->InsertNextPoint(sourceY);
    float sourceZ[3];
    vtkMath::Add(this->origin, this->zDirection, sourceZ);
    sourcePoints->InsertNextPoint(sourceZ);

    // Setup target points
    vtkSmartPointer<vtkPoints> targetPoints =
      vtkSmartPointer<vtkPoints>::New();

    targetPoints->InsertNextPoint(targetFrame.origin);
    float targetX[3];
    vtkMath::Add(targetFrame.origin, targetFrame.xDirection, targetX);
    targetPoints->InsertNextPoint(targetX);
    float targetY[3];
    vtkMath::Add(targetFrame.origin, targetFrame.yDirection, targetY);
    targetPoints->InsertNextPoint(targetY);
    float targetZ[3];
    vtkMath::Add(targetFrame.origin, targetFrame.zDirection, targetZ);
    targetPoints->InsertNextPoint(targetZ);

    landmarkTransform->SetSourceLandmarks(sourcePoints);
    landmarkTransform->SetTargetLandmarks(targetPoints);
    landmarkTransform->SetModeToRigidBody();
    landmarkTransform->Update();

    vtkMatrix4x4* M = landmarkTransform->GetMatrix();

    transform->SetMatrix(M);
  }

};

#endif