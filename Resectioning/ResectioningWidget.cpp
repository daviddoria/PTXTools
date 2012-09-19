/*=========================================================================
 *
 *  Copyright David Doria 2011 daviddoria@gmail.com
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

#include "ResectioningWidget.h"

// STL
#include <stdexcept>

// Eigen
#include <Eigen/Dense>

// ITK
#include "itkCastImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkVector.h"

// Qt
#include <QFileDialog>
#include <QIcon>
#include <QProgressDialog>
#include <QTextEdit>
#include <QtConcurrentRun>

// VTK
#include <vtkActor.h>
#include <vtkActor2D.h>
#include <vtkCamera.h>
#include <vtkCommand.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkFloatArray.h>
#include <vtkImageActor.h>
#include <vtkImageData.h>
#include <vtkImageSlice.h>
#include <vtkImageSliceMapper.h>
#include <vtkInteractorStyleImage.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkLookupTable.h>
#include <vtkMath.h>
#include <vtkPointData.h>
#include <vtkPointPicker.h>
#include <vtkProperty2D.h>
#include <vtkPolyDataMapper.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkStructuredGrid.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>

// Custom
#include "Custom3DStyle.h"
#include "CameraCalibration/CameraCalibration.h"
#include "PointSelectionStyle2D.h"
#include "PointSelectionStyle3D.h"
#include "Pane.h"
#include "Pane2D.h"
#include "Pane3D.h"
#include "PTXImage.h"
#include "PTXReader.h"
#include "Resectioning.h"
#include "ResectioningHelpers.h"
#include "Types.h"

// Submodules
#include "Helpers/Helpers.h"
#include "ITKHelpers/ITKHelpers.h"

void ResectioningWidget::on_actionHelp_activated()
{
  QTextEdit* help=new QTextEdit();
  
  help->setReadOnly(true);
  help->append("<h1>Image correspondences</h1>\
  Hold the right mouse button and drag to zoom in and out. <br/>\
  Hold the middle mouse button and drag to pan the image. <br/>\
  Click the left mouse button to select a keypoint.<br/> <p/>\
  <h1>Point cloud correspondences</h1>\
  Hold the left mouse button and drag to rotate the scene.<br/>\
  Hold the right mouse button and drag to zoom in and out. Hold the middle mouse button and drag to pan the scene. While holding control (CTRL), click the left mouse button to select a keypoint.<br/>\
  If you need to zoom in farther, hold shift while left clicking a point to change the camera's focal point to that point. You can reset the focal point by pressing 'r'.\
  <h1>Saving keypoints</h1>\
  The same number of keypoints must be selected in both the left and right panels before the points can be saved."
  );
  help->show();
}

void ResectioningWidget::on_actionQuit_activated()
{
  exit(0);
}

void ResectioningWidget::SharedConstructor()
{
  this->setupUi(this);

  ScannerLocationMarkerSource = vtkSmartPointer<vtkSphereSource>::New();
  ScannerLocationMarkerSource->Update();
  ScannerLocationMarkerMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  ScannerLocationMarkerMapper->SetInputData(ScannerLocationMarkerSource->GetOutput());
  ScannerLocationMarker = vtkSmartPointer<vtkActor>::New();
  ScannerLocationMarker->SetMapper(ScannerLocationMarkerMapper);
  
  CameraLocationMarkerSource = vtkSmartPointer<vtkSphereSource>::New();
  CameraLocationMarkerSource->Update();
  CameraLocationMarkerMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  CameraLocationMarkerMapper->SetInputData(CameraLocationMarkerSource->GetOutput());
  CameraLocationMarker = vtkSmartPointer<vtkActor>::New();
  CameraLocationMarker->SetMapper(CameraLocationMarkerMapper);

  Mesh = vtkSmartPointer<vtkPolyData>::New();
  InputMeshMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  InputMeshActor = vtkSmartPointer<vtkActor>::New();
  InputMeshActor->SetMapper(InputMeshMapper);
  InputMeshMapper->SetInputData(Mesh);

  OutputMeshMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  OutputMeshActor = vtkSmartPointer<vtkActor>::New();
  OutputMeshActor->SetMapper(OutputMeshMapper);
  OutputMeshMapper->SetInputData(Mesh);
  
  this->ProgressDialog = new QProgressDialog;
  this->ProgressDialog->setMinimum(0);
  this->ProgressDialog->setMaximum(0);
  this->ProgressDialog->setWindowModality(Qt::WindowModal);

  // Whenever the operation in another thread finishes, hide the progress bar.
  connect(&this->FutureWatcher, SIGNAL(finished()), this->ProgressDialog , SLOT(cancel()));

  this->PointCloudPane = new Pane3D(this->qvtkPointCloud);

  this->OutputPointCloudPane = new Pane3D(this->qvtkResult3D);
  
  this->ImagePane = new Pane2D(this->qvtkImage);

  this->ResultImagePane = new Pane2D(this->qvtkResultImage);

  this->ColorImage = FloatVectorImageType::New();
}

ResectioningWidget::ResectioningWidget(const std::string& imageFileName, const std::string& pointCloudFileName)
{
  SharedConstructor();

  LoadImage(imageFileName);

  LoadPTX(pointCloudFileName);
}

ResectioningWidget::ResectioningWidget(const std::string& imageFileName, const std::string& imageCorrespondenceFile,
                                       const std::string& pointCloudFileName,
                                       const std::string& pointCloudCorrespondenceFile)
{
  SharedConstructor();

  LoadImage(imageFileName);
  LoadCorrespondencesImage(imageCorrespondenceFile);

  LoadPTX(pointCloudFileName);
  LoadCorrespondencesPointCloud(pointCloudCorrespondenceFile);
}

ResectioningWidget::ResectioningWidget()
{
  SharedConstructor();
};

void ResectioningWidget::LoadCorrespondencesImage(const std::string& filename)
{
  std::string line;
  std::ifstream fin(filename.c_str());

  if(fin == NULL)
    {
    throw std::runtime_error("LoadCorrespondences2D Cannot open file.");
    }

  if(PointSelectionStyle2D::SafeDownCast(ImagePane->InteractorStyle))
    {
    PointSelectionStyle2D* selectionStyle = PointSelectionStyle2D::SafeDownCast(ImagePane->InteractorStyle);
    selectionStyle->RemoveAll();

    while(getline(fin, line))
      {
      std::stringstream ss;
      ss << line;
      double p[3];
      ss >> p[0] >> p[1];
      p[2] = 0;

      selectionStyle->AddNumber(p);
      }
    }
}

void ResectioningWidget::LoadCorrespondencesPointCloud(const std::string& filename)
{
  std::string line;
  std::ifstream fin(filename.c_str());

  if(fin == NULL)
    {
    throw std::runtime_error("LoadCorrespondences3D Cannot open file.");
    }

  if(PointSelectionStyle3D::SafeDownCast(PointCloudPane->InteractorStyle))
    {
    PointSelectionStyle3D* selectionStyle = PointSelectionStyle3D::SafeDownCast(PointCloudPane->InteractorStyle);

    selectionStyle->RemoveAll();

    while(getline(fin, line))
      {
      std::stringstream ss;
      ss << line;
      double p[3];
      ss >> p[0] >> p[1] >> p[2];

      selectionStyle->AddNumber(p);
      }
    }
}


void ResectioningWidget::SaveCorrespondencesImage(const std::string& filename)
{
  if(PointSelectionStyle2D::SafeDownCast(ImagePane->InteractorStyle))
    {
    PointSelectionStyle2D* selectionStyle = PointSelectionStyle2D::SafeDownCast(ImagePane->InteractorStyle);

    std::ofstream fout(filename.c_str());

    for(unsigned int i = 0; i < selectionStyle->GetNumberOfCorrespondences(); i++)
      {
      fout << selectionStyle->GetCorrespondence(i).x << " "
           << selectionStyle->GetCorrespondence(i).y << std::endl;

      }
    fout.close();
    }
}

void ResectioningWidget::SaveCorrespondencesPointCloud(const std::string& filename)
{
  if(PointSelectionStyle3D::SafeDownCast(PointCloudPane->InteractorStyle))
    {
    PointSelectionStyle3D* selectionStyle = PointSelectionStyle3D::SafeDownCast(PointCloudPane->InteractorStyle);

    std::ofstream fout(filename.c_str());

    for(unsigned int i = 0; i < selectionStyle->GetNumberOfCorrespondences(); i++)
      {
      fout << selectionStyle->GetCorrespondence(i).x << " "
          << selectionStyle->GetCorrespondence(i).y << " "
          << selectionStyle->GetCorrespondence(i).z << std::endl;
      }
    fout.close();
    }
}

void ResectioningWidget::LoadImage(const std::string& fileName)
{

/*
  QFileInfo fileInfo(fileName.toStdString().c_str());
  std::string extension = fileInfo.suffix().toStdString();
  std::cout << "extension: " << extension << std::endl;*/

  typedef itk::ImageFileReader<FloatVectorImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(fileName);
  reader->Update();

  ITKHelpers::DeepCopy(reader->GetOutput(), this->ColorImage.GetPointer());

  ImagePane->Image = reader->GetOutput();

  ResectioningHelpers::ITKImagetoVTKRGBImage(ImagePane->Image.GetPointer(), ImagePane->ImageData);

  ImagePane->ImageSliceMapper->SetInputData(ImagePane->ImageData);
  ImagePane->ImageSlice->SetMapper(ImagePane->ImageSliceMapper);
  
  // Add Actor to renderer
  //pane->Renderer->AddActor(pane->ImageActor);
  ImagePane->Renderer->AddActor(ImagePane->ImageSlice);
  ImagePane->Renderer->ResetCamera();

  vtkSmartPointer<vtkPointPicker> pointPicker = vtkSmartPointer<vtkPointPicker>::New();
  ImagePane->qvtkWidget->GetRenderWindow()->GetInteractor()->SetPicker(pointPicker);

  ImagePane->InteractorStyle = PointSelectionStyle2D::New();
  ImagePane->InteractorStyle->SetCurrentRenderer(ImagePane->Renderer);

  //pane->SelectionStyle->Initialize();
  ImagePane->qvtkWidget->GetRenderWindow()->GetInteractor()->
            SetInteractorStyle(static_cast<PointSelectionStyle2D*>(ImagePane->InteractorStyle));

  ImagePane->Renderer->ResetCamera();

  ImagePane->qvtkWidget->GetRenderWindow()->Render();
}

void ResectioningWidget::on_action_Image_FlipHorizontally_activated()
{
  ImagePane->FlipHorizontally();
}

void ResectioningWidget::on_action_Image_FlipVertically_activated()
{
  ImagePane->FlipVertically();
}

void ResectioningWidget::LoadPTX(const std::string& fileName)
{
  PTXReader reader;
  reader.SetFileName(fileName.c_str());

  // Read in the same thread
//   reader.Read();
//   PTXImage ptxImage = reader.GetOutput();

  // Read the file in a different thread because it could take a while.
  QFuture<void> readerFuture = QtConcurrent::run(&reader, &PTXReader::Read);
  this->FutureWatcher.setFuture(readerFuture);
  this->ProgressDialog->setLabelText("Opening PTX file...");
  this->ProgressDialog->exec();

  // Since we have displayed a modal dialog, it will wait here without needing to call WaitForFinished

  this->PTX = reader.GetOutput();
  this->OriginalPTX = reader.GetOutput();

  vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
  this->PTX.CreatePointCloud(polyData);

  vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();
  this->PTX.CreateStructuredGrid(structuredGrid);
  ResectioningHelpers::StructuredGridToPolyData(structuredGrid, Mesh);

  vtkSmartPointer<vtkLookupTable> lookupTable = vtkSmartPointer<vtkLookupTable>::New();
  //lookupTable->SetTableRange(0.0, 10.0);
  //lookupTable->SetHueRange(0, .5);
  //lookupTable->SetHueRange(.5, 1);
  lookupTable->SetHueRange(0, 1);

  // If the cloud has an "Intensity" array, use it for the initial coloring
  vtkFloatArray* intensityArray = vtkFloatArray::SafeDownCast(
                                  polyData->GetPointData()->GetArray("Intensity"));
  if(intensityArray)
    {
    polyData->GetPointData()->SetActiveScalars("Intensity");

    float range[2];
    intensityArray->GetValueRange(range);

    lookupTable->SetTableRange(range[0], range[1]);
    }

  PointCloudPane->PointCloudMapper->SetInputData(polyData);
  PointCloudPane->PointCloudMapper->SetLookupTable(lookupTable);
  PointCloudPane->PointCloudActor->SetMapper(PointCloudPane->PointCloudMapper);
  PointCloudPane->PointCloudActor->GetProperty()->SetRepresentationToPoints();

  // Add Actor to renderer
  PointCloudPane->Renderer->AddActor(PointCloudPane->PointCloudActor);
  PointCloudPane->Renderer->ResetCamera();

  vtkSmartPointer<vtkPointPicker> pointPicker = vtkSmartPointer<vtkPointPicker>::New();

  pointPicker->PickFromListOn();
  pointPicker->AddPickList(PointCloudPane->PointCloudActor);
  PointCloudPane->qvtkWidget->GetRenderWindow()->GetInteractor()->SetPicker(pointPicker);

  PointCloudPane->InteractorStyle = PointSelectionStyle3D::New();
  PointCloudPane->InteractorStyle->SetCurrentRenderer(PointCloudPane->Renderer);
  dynamic_cast<PointSelectionStyle3D*>(PointCloudPane->InteractorStyle)->Initialize();
  
  static_cast<PointSelectionStyle3D*>(PointCloudPane->InteractorStyle)->Data = polyData;
  PointCloudPane->qvtkWidget->GetRenderWindow()->GetInteractor()->SetInteractorStyle(
        static_cast<PointSelectionStyle3D*>(PointCloudPane->InteractorStyle));

  PointCloudPane->Renderer->ResetCamera();

  std::cout << "Computing average spacing..." << std::endl;
  //float averageSpacing = ResectioningHelpers::ComputeAverageSpacing(polyData->GetPoints(), 100000);
  unsigned int pointsToUse = 1000;
  QFuture<float> spacingFuture = QtConcurrent::run(ResectioningHelpers::ComputeAverageSpacing,
                                                  polyData->GetPoints(), pointsToUse);
  this->FutureWatcher.setFuture(spacingFuture);
  this->ProgressDialog->setLabelText("Computing average spacing...");
  this->ProgressDialog->exec();

  float averageSpacing = *spacingFuture.begin();
  std::cout << "Done computing average spacing: " << averageSpacing << std::endl;

  dynamic_cast<PointSelectionStyle3D*>(PointCloudPane->InteractorStyle)->SetMarkerRadius(averageSpacing * 10.0);
}

void ResectioningWidget::LoadVTP(const std::string& fileName)
{

/*
  QFileInfo fileInfo(fileName.toStdString().c_str());
  std::string extension = fileInfo.suffix().toStdString();
  std::cout << "extension: " << extension << std::endl;
  */

  vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
  reader->SetFileName(fileName.c_str());

  //reader->Update();
  // Start the computation.
  QFuture<void> future = QtConcurrent::run(reader.GetPointer(),
                            static_cast<void(vtkXMLPolyDataReader::*)()>(&vtkXMLPolyDataReader::Update));
  this->FutureWatcher.setFuture(future);
  this->ProgressDialog->setMinimum(0);
  this->ProgressDialog->setMaximum(0);
  this->ProgressDialog->setLabelText("Opening file...");
  this->ProgressDialog->setWindowModality(Qt::WindowModal);
  this->ProgressDialog->exec();

  vtkSmartPointer<vtkLookupTable> lookupTable = vtkSmartPointer<vtkLookupTable>::New();
  //lookupTable->SetTableRange(0.0, 10.0);
  //lookupTable->SetHueRange(0, .5);
  //lookupTable->SetHueRange(.5, 1);
  lookupTable->SetHueRange(0, 1);

  // If the cloud has an "Intensity" array, use it for the initial coloring
  vtkFloatArray* intensityArray = vtkFloatArray::SafeDownCast(
                                  reader->GetOutput()->GetPointData()->GetArray("Intensity"));
  if(intensityArray)
    {
    reader->GetOutput()->GetPointData()->SetActiveScalars("Intensity");

    float range[2];
    intensityArray->GetValueRange(range);

    lookupTable->SetTableRange(range[0], range[1]);
    }

  PointCloudPane->PointCloudMapper->SetInputConnection(reader->GetOutputPort());
  PointCloudPane->PointCloudMapper->SetLookupTable(lookupTable);
  PointCloudPane->PointCloudActor->SetMapper(PointCloudPane->PointCloudMapper);
  PointCloudPane->PointCloudActor->GetProperty()->SetRepresentationToPoints();

  // Add Actor to renderer
  PointCloudPane->Renderer->AddActor(PointCloudPane->PointCloudActor);
  PointCloudPane->Renderer->ResetCamera();

  vtkSmartPointer<vtkPointPicker> pointPicker = vtkSmartPointer<vtkPointPicker>::New();

  pointPicker->PickFromListOn();
  pointPicker->AddPickList(PointCloudPane->PointCloudActor);
  PointCloudPane->qvtkWidget->GetRenderWindow()->GetInteractor()->SetPicker(pointPicker);
  PointCloudPane->InteractorStyle = PointSelectionStyle3D::New();
  PointCloudPane->InteractorStyle->SetCurrentRenderer(PointCloudPane->Renderer);
  dynamic_cast<PointSelectionStyle3D*>(PointCloudPane->InteractorStyle)->Initialize();
  dynamic_cast<PointSelectionStyle3D*>(PointCloudPane->InteractorStyle)->Data = reader->GetOutput();
  PointCloudPane->qvtkWidget->GetRenderWindow()->GetInteractor()->SetInteractorStyle(
         dynamic_cast<PointSelectionStyle3D*>(PointCloudPane->InteractorStyle));

  PointCloudPane->Renderer->ResetCamera();

  std::cout << "Computing average spacing..." << std::endl;
  float averageSpacing = ResectioningHelpers::ComputeAverageSpacing(reader->GetOutput()->GetPoints(), 100000);
  std::cout << "Done computing average spacing: " << averageSpacing << std::endl;

  dynamic_cast<PointSelectionStyle3D*>(PointCloudPane->InteractorStyle)->SetMarkerRadius(averageSpacing * 10.0);
}

void ResectioningWidget::on_action_Image_LoadCorrespondences_activated()
{
  // Get a filename to open
  QString fileName = QFileDialog::getOpenFileName(this, "Open File", ".",
                                                  "Image Files (*.jpg *.jpeg *.bmp *.png *.mha)");

  std::cout << "Got filename: " << fileName.toStdString() << std::endl;
  if(fileName.toStdString().empty())
    {
    std::cout << "Filename was empty." << std::endl;
    return;
    }

  LoadCorrespondencesImage(fileName.toStdString());
}

void ResectioningWidget::on_action_PointCloud_LoadCorrespondences_activated()
{
  // Get a filename to open
  QString fileName = QFileDialog::getOpenFileName(this, "Open File", ".",
                                                  "Image Files (*.jpg *.jpeg *.bmp *.png *.mha)");

  std::cout << "Got filename: " << fileName.toStdString() << std::endl;
  if(fileName.toStdString().empty())
    {
    std::cout << "Filename was empty." << std::endl;
    return;
    }
  LoadCorrespondencesPointCloud(fileName.toStdString());
}

void ResectioningWidget::on_action_Image_Open_activated()
{
  // Get a filename to open
  QString fileName = QFileDialog::getOpenFileName(this, "Open File", ".",
                                                  "Image Files (*.jpg *.jpeg *.bmp *.png *.mha)");

  std::cout << "Got filename: " << fileName.toStdString() << std::endl;
  if(fileName.toStdString().empty())
    {
    std::cout << "Filename was empty." << std::endl;
    return;
    }

  LoadImage(fileName.toStdString());
}

void ResectioningWidget::on_action_PointCloud_OpenVTP_activated()
{
  // Get a filename to open
  QString fileName = QFileDialog::getOpenFileName(this, "Open Point Cloud", ".", "Files (*.vtp)");

  std::cout << "Got filename: " << fileName.toStdString() << std::endl;
  if(fileName.toStdString().empty())
    {
    std::cout << "Filename was empty." << std::endl;
    return;
    }

  LoadVTP(fileName.toStdString());
  std::cout << "Done loading point cloud." << std::endl;
}

void ResectioningWidget::on_action_PointCloud_OpenPTX_activated()
{
  // Get a filename to open
  QString fileName = QFileDialog::getOpenFileName(this, "Open Point Cloud", ".", "Files (*.vtp)");

  std::cout << "Got filename: " << fileName.toStdString() << std::endl;
  if(fileName.toStdString().empty())
    {
    std::cout << "Filename was empty." << std::endl;
    return;
    }

  LoadPTX(fileName.toStdString());
  std::cout << "Done loading point cloud." << std::endl;
}


void ResectioningWidget::on_action_Image_SaveCorrespondences_activated()
{
  QString fileName = QFileDialog::getSaveFileName(this, "Save File", ".", "Text Files (*.txt)");
  std::cout << "Got filename: " << fileName.toStdString() << std::endl;
  if(fileName.toStdString().empty())
    {
    std::cout << "Filename was empty." << std::endl;
    return;
    }
  SaveCorrespondencesImage(fileName.toStdString());
  
}

void ResectioningWidget::on_action_PointCloud_SaveCorrespondences_activated()
{
  QString fileName = QFileDialog::getSaveFileName(this, "Save File", ".", "Text Files (*.txt)");
  std::cout << "Got filename: " << fileName.toStdString() << std::endl;
  if(fileName.toStdString().empty())
    {
    std::cout << "Filename was empty." << std::endl;
    return;
    }
  SaveCorrespondencesPointCloud(fileName.toStdString());
}

void ResectioningWidget::on_btnDeleteLastCorrespondencePointCloud_clicked()
{
  dynamic_cast<PointSelectionStyle3D*>(this->PointCloudPane->InteractorStyle)->DeleteLastCorrespondence();
  this->PointCloudPane->Refresh();
}

void ResectioningWidget::on_btnDeleteAllCorrespondencesPointCloud_clicked()
{
  dynamic_cast<PointSelectionStyle3D*>(this->PointCloudPane->InteractorStyle)->RemoveAll();
  this->PointCloudPane->Refresh();
  //static_cast<PointSelectionStyle2D*>(pane->SelectionStyle)->RemoveAll();
  //this->qvtkWidgetLeft->GetRenderWindow()->Render();
}

void ResectioningWidget::on_btnDeleteLastCorrespondenceImage_clicked()
{
  dynamic_cast<PointSelectionStyle2D*>(this->ImagePane->InteractorStyle)->DeleteLastCorrespondence();
  this->ImagePane->Refresh();
}

void ResectioningWidget::on_btnDeleteAllCorrespondencesImage_clicked()
{
  dynamic_cast<PointSelectionStyle2D*>(this->ImagePane->InteractorStyle)->RemoveAll();
  this->ImagePane->Refresh();
  //static_cast<PointSelectionStyle2D*>(pane->SelectionStyle)->RemoveAll();
  //this->qvtkWidgetRight->GetRenderWindow()->Render();
}

Eigen::MatrixXd ResectioningWidget::ComputeP()
{
  CameraCalibration::Point2DVector points2D;
  std::cout << "There are " << dynamic_cast<PointSelectionStyle2D*>(ImagePane->InteractorStyle)->GetNumberOfCorrespondences()
            << " correspondences in the image pane." << std::endl;

  for(vtkIdType pointId = 0;
      pointId < static_cast<vtkIdType>(dynamic_cast<PointSelectionStyle2D*>(ImagePane->InteractorStyle)->GetNumberOfCorrespondences()); ++pointId)
  {
    Coord3D coord = dynamic_cast<PointSelectionStyle2D*>(ImagePane->InteractorStyle)->GetCorrespondence(pointId);
    points2D.push_back(Eigen::Vector2d (coord.x, coord.y));
  }

  CameraCalibration::Point3DVector points3D;
  std::cout << "There are " << dynamic_cast<PointSelectionStyle3D*>(PointCloudPane->InteractorStyle)->GetNumberOfCorrespondences()
            << " correspondences in the point cloud pane." << std::endl;
  for(vtkIdType pointId = 0;
      pointId < static_cast<vtkIdType>(dynamic_cast<PointSelectionStyle3D*>(PointCloudPane->InteractorStyle)->GetNumberOfCorrespondences()); ++pointId)
  {
    Coord3D coord = dynamic_cast<PointSelectionStyle3D*>(PointCloudPane->InteractorStyle)->GetCorrespondence(pointId);
    points3D.push_back(Eigen::Vector3d (coord.x, coord.y, coord.z));
  }

  std::cout << "There are " << points2D.size() << " 2D points and "
            << points3D.size() << " 3D points to use with DLT." << std::endl;
  Eigen::MatrixXd P = CameraCalibration::ComputeP_NormalizedDLT(points2D, points3D);

  return P;
}

void ResectioningWidget::on_btnComputeP_clicked()
{
  Eigen::MatrixXd P = ComputeP();

  //double cameraLocation[3] = {P(0, 3), P(1,3), P(2,3)};

  Eigen::VectorXd C = CameraCalibration::GetCameraCenter(P);

  double cameraLocation[3] = {C[0], C[1], C[2]};

  std::cout << "cameraLocation: " << cameraLocation[0] << " " << cameraLocation[1] << " " << cameraLocation[2] << std::endl;
  CameraLocationMarker->SetPosition(cameraLocation);

  PointCloudPane->Renderer->AddActor(CameraLocationMarker);

  PointCloudPane->Renderer->AddActor(ScannerLocationMarker);
  
  PointCloudPane->qvtkWidget->GetRenderWindow()->Render();
}

void ResectioningWidget::on_btnResection_clicked()
{
  Eigen::MatrixXd P = ComputeP();

  typedef itk::Image<itk::RGBPixel<unsigned char>, 2> RGBImageType;
  RGBImageType::Pointer rgbImage = RGBImageType::New();

  //Helpers::WriteImage(this->ColorImage.GetPointer(), "colorImage.mha");

  ITKHelpers::VectorImageToRGBImage(this->ColorImage.GetPointer(), rgbImage.GetPointer());

  //Helpers::WriteImage(rgbImage.GetPointer(), "rgbImage.png");

  // Compute the resectioning in the same thread
  //this->PTX = Resectioning::ResectionSmart(P, this->OriginalPTX, rgbImage.GetPointer());
  //this->PTX = Resectioning::ResectionNaive(P, this->OriginalPTX, rgbImage.GetPointer());

  // Compute the resectioning in a different thread
//   QFuture<PTXImage> resectionFuture = QtConcurrent::run(Resectioning::Resection_ProjectionSorting, P,
//                                                         this->OriginalPTX, rgbImage.GetPointer());

  QFuture<PTXImage> resectionFuture = QtConcurrent::run(Resectioning::Resection_ProjectionSorting, P,
                                                        this->OriginalPTX, rgbImage.GetPointer());

  this->FutureWatcher.setFuture(resectionFuture);
  this->ProgressDialog->setLabelText("Resectioning...");
  this->ProgressDialog->exec();
  this->PTX = *resectionFuture.begin();

  ShowResultImage();
  ShowResultPointCloud();
}

void ResectioningWidget::on_btnResectionMesh_clicked()
{
  // this->PTX.GetMesh(this->Mesh);

  Eigen::MatrixXd P = ComputeP();

  typedef itk::Image<itk::RGBPixel<unsigned char>, 2> RGBImageType;
  RGBImageType::Pointer rgbImage = RGBImageType::New();

  //Helpers::WriteImage(this->ColorImage.GetPointer(), "colorImage.mha");
  
  ITKHelpers::VectorImageToRGBImage(this->ColorImage.GetPointer(), rgbImage.GetPointer());

  //Helpers::WriteImage(rgbImage.GetPointer(), "rgbImage.png");
  
  // Compute the resectioning in the same thread
  //this->PTX = Resectioning::ResectionSmart(P, this->OriginalPTX, rgbImage.GetPointer());
  //this->PTX = Resectioning::ResectionNaive(P, this->OriginalPTX, rgbImage.GetPointer());

  // Compute the resectioning in a different thread
//   QFuture<PTXImage> resectionFuture = QtConcurrent::run(Resectioning::Resection_ProjectionSorting, P,
//                                                         this->OriginalPTX, rgbImage.GetPointer());

  if(PTX.GetMesh()->GetNumberOfPoints() == 0)
    {
    float maxEdgeLength = 1.0f;
    QFuture<void> triangulationFuture = QtConcurrent::run(&PTX, &PTXImage::ComputeMesh, maxEdgeLength);
    this->FutureWatcher.setFuture(triangulationFuture);
    this->ProgressDialog->setLabelText("Triangulating...");
    this->ProgressDialog->exec();
    }

  QFuture<PTXImage> resectionFuture = QtConcurrent::run(Resectioning::Resection_MeshIntersection, P,
                                                        this->OriginalPTX, rgbImage.GetPointer());

  this->FutureWatcher.setFuture(resectionFuture);
  this->ProgressDialog->setLabelText("Resectioning using mesh intersection...");
  this->ProgressDialog->exec();
  this->PTX = *resectionFuture.begin();

  ShowResultImage();
}

void ResectioningWidget::ShowResultPointCloud()
{
  vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
  this->PTX.CreatePointCloud(polyData);

  this->OutputPointCloudPane->SetPolyData(polyData);

  this->OutputPointCloudPane->Renderer->ResetCamera();

  OutputPointCloudPane->InteractorStyle = Custom3DStyle::New();
  OutputPointCloudPane->InteractorStyle->SetCurrentRenderer(OutputPointCloudPane->Renderer);

  OutputPointCloudPane->qvtkWidget->GetRenderWindow()->GetInteractor()->SetInteractorStyle(OutputPointCloudPane->InteractorStyle);

  qvtkResult3D->GetRenderWindow()->Render();
}

void ResectioningWidget::ShowResultImage()
{
  PTXImage::RGBImageType::Pointer rgbimage = PTXImage::RGBImageType::New();
  this->PTX.CreateRGBImage(rgbimage.GetPointer());

  typedef itk::VectorImage<float, 2> VectorImageType;
  VectorImageType::Pointer image = VectorImageType::New();

  ITKHelpers::RGBImageToVectorImage(rgbimage.GetPointer(), image.GetPointer());

  ResultImagePane->Image = image;

  ResectioningHelpers::ITKImagetoVTKRGBImage(ResultImagePane->Image.GetPointer(), ResultImagePane->ImageData);

  ResultImagePane->ImageSliceMapper->SetInputData(ResultImagePane->ImageData);
  ResultImagePane->ImageSlice->SetMapper(ResultImagePane->ImageSliceMapper);

  // Add Actor to renderer
  ResultImagePane->Renderer->AddActor(ResultImagePane->ImageSlice);

  ResultImagePane->InteractorStyle = vtkInteractorStyleImage::New();
  ResultImagePane->InteractorStyle->SetCurrentRenderer(ResultImagePane->Renderer);
  ResultImagePane->qvtkWidget->GetRenderWindow()->GetInteractor()->SetInteractorStyle(ResultImagePane->InteractorStyle);

  ResultImagePane->Renderer->ResetCamera();

  qvtkResultImage->GetRenderWindow()->Render();
}

void ResectioningWidget::on_action_Export_ResultPTX_activated()
{
  QString fileName = QFileDialog::getSaveFileName(this, "Save PTX", ".",
                                                  "PTX Files (*.ptx)");

  std::cout << "Got filename: " << fileName.toStdString() << std::endl;
  if(fileName.toStdString().empty())
    {
    std::cout << "Filename was empty." << std::endl;
    return;
    }

  QFuture<void> writerFuture = QtConcurrent::run(&PTX, &PTXImage::WritePTX, fileName.toStdString());
  this->FutureWatcher.setFuture(writerFuture);
  this->ProgressDialog->setLabelText("Writing PTX file...");
  this->ProgressDialog->exec();

  //this->PTX.WritePTX(fileName.toStdString());
}

void ResectioningWidget::on_action_Export_ResultRGB_activated()
{
  QString fileName = QFileDialog::getSaveFileName(this, "Save Image", ".",
                                                  "Image Files (*.jpg *.jpeg *.bmp *.png *.mha)");

  std::cout << "Got filename: " << fileName.toStdString() << std::endl;
  if(fileName.toStdString().empty())
    {
    std::cout << "Filename was empty." << std::endl;
    return;
    }
  this->PTX.WriteRGBImage(fileName.toStdString());
}

void ResectioningWidget::on_btnWriteMesh_clicked()
{
  QString fileName = QFileDialog::getSaveFileName(this, "Save VTP", ".",
                                                  "VTP Files (*.vtp)");
  std::cout << "Got filename: " << fileName.toStdString() << std::endl;
  if(fileName.toStdString().empty())
    {
    std::cout << "Filename was empty." << std::endl;
    return;
    }
    
  if(PTX.GetMesh()->GetNumberOfPoints() == 0)
    {
    QFuture<void> triangulationFuture = QtConcurrent::run(&PTX, &PTXImage::ComputeMesh, 1.0f);
    this->FutureWatcher.setFuture(triangulationFuture);
    this->ProgressDialog->setLabelText("Triangulating...");
    this->ProgressDialog->exec();
    }

  vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName(fileName.toStdString().c_str());
  writer->SetInputData(PTX.GetMesh());
  writer->Write();
}

void ResectioningWidget::on_chkInputMesh_clicked()
{
  if(chkInputMesh->isChecked())
  {
    PointCloudPane->Renderer->AddActor(InputMeshActor);
  }
  else
  {
    PointCloudPane->Renderer->RemoveActor(InputMeshActor);
  }

  PointCloudPane->Refresh();
}

void ResectioningWidget::on_chkOutputMesh_clicked()
{
  if(chkOutputMesh->isChecked())
  {
    OutputPointCloudPane->Renderer->AddActor(OutputMeshActor);
  }
  else
  {
    OutputPointCloudPane->Renderer->RemoveActor(OutputMeshActor);
  }
  OutputPointCloudPane->Refresh();
}
