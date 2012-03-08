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

#include "ui_SelectCorrespondencesWidget.h"
#include "ResectioningWidget.h"

// STL
#include <stdexcept>

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
#include <vtkInteractorStyleImage.h>
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
#include <vtkImageSliceMapper.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkXMLPolyDataReader.h>

// Custom
#include "Helpers.h"
#include "ResectioningHelpers.h"
#include "Types.h"
#include "PointSelectionStyle2D.h"
#include "PointSelectionStyle3D.h"
#include "Pane.h"
#include "Pane2D.h"
#include "Pane3D.h"

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
  LeftPane = NULL;
  RightPane = NULL;
  ProgressDialog = new QProgressDialog;

  this->setupUi(this);

  connect(&this->FutureWatcher, SIGNAL(finished()), this->ProgressDialog , SLOT(cancel()));

  // Setup icons
  QIcon openIcon = QIcon::fromTheme("document-open");
  QIcon saveIcon = QIcon::fromTheme("document-save");

  // Setup left toolbar
  actionOpenImageLeft->setIcon(openIcon);
  this->toolBar_left->addAction(actionOpenImageLeft);

  actionOpenPointCloudLeft->setIcon(openIcon);
  this->toolBar_left->addAction(actionOpenPointCloudLeft);

  actionSavePointsLeft->setIcon(saveIcon);
  this->toolBar_left->addAction(actionSavePointsLeft);

  actionLoadPointsLeft->setIcon(openIcon);
  this->toolBar_left->addAction(actionLoadPointsLeft);

  // Setup right toolbar
  actionOpenImageRight->setIcon(openIcon);
  this->toolBar_right->addAction(actionOpenImageRight);

  actionOpenPointCloudRight->setIcon(openIcon);
  this->toolBar_right->addAction(actionOpenPointCloudRight);

  actionSavePointsRight->setIcon(saveIcon);
  this->toolBar_right->addAction(actionSavePointsRight);

  actionLoadPointsRight->setIcon(openIcon);
  this->toolBar_right->addAction(actionLoadPointsRight);
}

ResectioningWidget::ResectioningWidget(const std::string& imageFileName, const std::string& pointCloudFileName)
{
  SharedConstructor();

  this->LeftPane = new Pane2D(this->qvtkWidgetLeft);
  LoadImage(LeftPane, imageFileName);

  this->RightPane = new Pane3D(this->qvtkWidgetRight);
  LoadPointCloud(RightPane, pointCloudFileName);
}

// Constructor
ResectioningWidget::ResectioningWidget()
{
  SharedConstructor();

};

void ResectioningWidget::LoadPoints(Pane* const pane, const std::string& fileName)
{
  // This function reads existing correspondences from a plain text file and displays them in the left pane.

  if(!this->LeftPane)
    {
    std::cerr << "Cannot load points unless an image or point cloud is loaded!" << std::endl;
    return;
    }

  if(dynamic_cast<Pane2D*>(pane))
    {
    LoadPoints2D(dynamic_cast<Pane2D*>(pane), fileName);
    }
  else if(dynamic_cast<Pane3D*>(pane))
    {
    LoadPoints3D(dynamic_cast<Pane3D*>(pane), fileName);
    }
}


void ResectioningWidget::LoadPoints2D(Pane2D* const pane, const std::string& filename)
{
  std::string line;
  std::ifstream fin(filename.c_str());

  if(fin == NULL)
    {
    throw std::runtime_error("LoadPoints2D Cannot open file.");
    }

  pane->SelectionStyle->RemoveAll();

  while(getline(fin, line))
    {
    std::stringstream ss;
    ss << line;
    double p[3];
    ss >> p[0] >> p[1];
    p[2] = 0;

    pane->SelectionStyle->AddNumber(p);
    }
}

void ResectioningWidget::LoadPoints3D(Pane3D* const pane, const std::string& filename)
{
  std::string line;
  std::ifstream fin(filename.c_str());

  if(fin == NULL)
    {
    throw std::runtime_error("LoadPoints3D Cannot open file.");
    }

  pane->SelectionStyle->RemoveAll();

  while(getline(fin, line))
    {
    std::stringstream ss;
    ss << line;
    double p[3];
    ss >> p[0] >> p[1] >> p[2];

    pane->SelectionStyle->AddNumber(p);
    }
}

void ResectioningWidget::SavePoints(Pane* const pane)
{
  if(!pane)
    {
    std::cerr << "You must have loaded and selected points from this pane!" << std::endl;
    return;
    }

  QString fileName = QFileDialog::getSaveFileName(this, "Save File", ".", "Text Files (*.txt)");
  std::cout << "Got filename: " << fileName.toStdString() << std::endl;
  if(fileName.toStdString().empty())
    {
    std::cout << "Filename was empty." << std::endl;
    return;
    }

  if(dynamic_cast<Pane2D*>(pane))
    {
    SavePoints2D(dynamic_cast<Pane2D*>(pane), fileName.toStdString());
    }
  else if(dynamic_cast<Pane3D*>(pane))
    {
    SavePoints3D(dynamic_cast<Pane3D*>(pane), fileName.toStdString());
    }
  
}

void ResectioningWidget::SavePoints2D(Pane2D* const pane, const std::string& filename)
{
  std::ofstream fout(filename.c_str());

  for(unsigned int i = 0; i < pane->SelectionStyle->GetNumberOfCorrespondences(); i++)
    {
    fout << pane->SelectionStyle->GetCorrespondence(i).x << " "
         << pane->SelectionStyle->GetCorrespondence(i).y << std::endl;

    }
  fout.close();
}

void ResectioningWidget::SavePoints3D(Pane3D* const pane, const std::string& filename)
{
  std::ofstream fout(filename.c_str());

  for(unsigned int i = 0; i < pane->SelectionStyle->GetNumberOfCorrespondences(); i++)
    {
    fout << pane->SelectionStyle->GetCorrespondence(i).x << " "
         << pane->SelectionStyle->GetCorrespondence(i).y << " "
         << pane->SelectionStyle->GetCorrespondence(i).z << std::endl;
    }
  fout.close();
}

void ResectioningWidget::LoadImage(Pane* const inputPane, const std::string& fileName)
{

/*
  QFileInfo fileInfo(fileName.toStdString().c_str());
  std::string extension = fileInfo.suffix().toStdString();
  std::cout << "extension: " << extension << std::endl;*/
  
  Pane2D* pane = static_cast<Pane2D*>(inputPane);

  if(!pane)
  {
    throw std::runtime_error("LoadImage: inputPane is NULL!");
  }
  typedef itk::ImageFileReader<FloatVectorImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(fileName);
  reader->Update();

  pane->Image = reader->GetOutput();

  if(this->chkRGB->isChecked())
    {
    ResectioningHelpers::ITKImagetoVTKRGBImage(pane->Image.GetPointer(), pane->ImageData);
    }
  else
    {
    ResectioningHelpers::ITKImagetoVTKMagnitudeImage(pane->Image.GetPointer(), pane->ImageData);
    }

  pane->ImageSliceMapper->SetInputConnection(pane->ImageData->GetProducerPort());
  pane->ImageSlice->SetMapper(pane->ImageSliceMapper);
  
  // Add Actor to renderer
  //pane->Renderer->AddActor(pane->ImageActor);
  pane->Renderer->AddActor(pane->ImageSlice);
  pane->Renderer->ResetCamera();

  vtkSmartPointer<vtkPointPicker> pointPicker = vtkSmartPointer<vtkPointPicker>::New();
  pane->qvtkWidget->GetRenderWindow()->GetInteractor()->SetPicker(pointPicker);
  pane->SelectionStyle = PointSelectionStyle2D::New();
  pane->SelectionStyle->SetCurrentRenderer(pane->Renderer);
  //pane->SelectionStyle->Initialize();
  pane->qvtkWidget->GetRenderWindow()->GetInteractor()->SetInteractorStyle(static_cast<PointSelectionStyle2D*>(pane->SelectionStyle));

  pane->Renderer->ResetCamera();

  pane->qvtkWidget->GetRenderWindow()->Render();
}

void ResectioningWidget::on_actionFlipLeftHorizontally_activated()
{
  if(dynamic_cast<Pane2D*>(this->LeftPane))
    {
    static_cast<Pane2D*>(this->LeftPane)->FlipHorizontally();
    }
  else
    {
    std::cerr << "Cannot flip a point cloud!" << std::endl;
    }
}

void ResectioningWidget::on_actionFlipLeftVertically_activated()
{
  if(dynamic_cast<Pane2D*>(this->LeftPane))
    {
    static_cast<Pane2D*>(this->LeftPane)->FlipVertically();
    }
  else
    {
    std::cerr << "Cannot flip a point cloud!" << std::endl;
    }
}

void ResectioningWidget::on_actionFlipRightHorizontally_activated()
{
  if(dynamic_cast<Pane2D*>(this->RightPane))
    {
    static_cast<Pane2D*>(this->RightPane)->FlipHorizontally();
    }
  else
    {
    std::cerr << "Cannot flip a point cloud!" << std::endl;
    }
}

void ResectioningWidget::on_actionFlipRightVertically_activated()
{
  if(dynamic_cast<Pane2D*>(this->RightPane))
    {
    static_cast<Pane2D*>(this->RightPane)->FlipVertically();
    }
  else
    {
    std::cerr << "Cannot flip a point cloud!" << std::endl;
    }
}

void ResectioningWidget::LoadPointCloud(Pane* const inputPane, const std::string& fileName)
{

/*
  QFileInfo fileInfo(fileName.toStdString().c_str());
  std::string extension = fileInfo.suffix().toStdString();
  std::cout << "extension: " << extension << std::endl;
  */
  Pane3D* pane = static_cast<Pane3D*>(inputPane);

  if(!pane)
  {
    throw std::runtime_error("LoadImage: inputPane is NULL!");
  }
  
  vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
  reader->SetFileName(fileName.c_str());

  //reader->Update();
  // Start the computation.
  QFuture<void> future = QtConcurrent::run(reader.GetPointer(), static_cast<void(vtkXMLPolyDataReader::*)()>(&vtkXMLPolyDataReader::Update));
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
  vtkFloatArray* intensityArray = vtkFloatArray::SafeDownCast(reader->GetOutput()->GetPointData()->GetArray("Intensity"));
  if(intensityArray)
    {
    reader->GetOutput()->GetPointData()->SetActiveScalars("Intensity");

    float range[2];
    intensityArray->GetValueRange(range);

    lookupTable->SetTableRange(range[0], range[1]);
    }

  pane->PointCloudMapper->SetInputConnection(reader->GetOutputPort());
  pane->PointCloudMapper->SetLookupTable(lookupTable);
  pane->PointCloudActor->SetMapper(pane->PointCloudMapper);
  pane->PointCloudActor->GetProperty()->SetRepresentationToPoints();

  // Add Actor to renderer
  pane->Renderer->AddActor(pane->PointCloudActor);
  pane->Renderer->ResetCamera();

  vtkSmartPointer<vtkPointPicker> pointPicker = vtkSmartPointer<vtkPointPicker>::New();

  pointPicker->PickFromListOn();
  pointPicker->AddPickList(pane->PointCloudActor);
  pane->qvtkWidget->GetRenderWindow()->GetInteractor()->SetPicker(pointPicker);
  pane->SelectionStyle = PointSelectionStyle3D::New();
  pane->SelectionStyle->SetCurrentRenderer(pane->Renderer);
  pane->SelectionStyle->Initialize();
  static_cast<PointSelectionStyle3D*>(pane->SelectionStyle)->Data = reader->GetOutput();
  pane->qvtkWidget->GetRenderWindow()->GetInteractor()->SetInteractorStyle(static_cast<PointSelectionStyle3D*>(pane->SelectionStyle));

  pane->Renderer->ResetCamera();

  std::cout << "Computing average spacing..." << std::endl;
  float averageSpacing = ResectioningHelpers::ComputeAverageSpacing(reader->GetOutput()->GetPoints(), 100000);
  std::cout << "Done computing average spacing: " << averageSpacing << std::endl;

  static_cast<PointSelectionStyle3D*>(pane->SelectionStyle)->SetMarkerRadius(averageSpacing * 10.0);
}

void ResectioningWidget::on_actionLoadPointsLeft_activated()
{
  // Get a filename to open
  QString fileName = QFileDialog::getOpenFileName(this, "Open File", ".", "Image Files (*.jpg *.jpeg *.bmp *.png *.mha)");

  std::cout << "Got filename: " << fileName.toStdString() << std::endl;
  if(fileName.toStdString().empty())
    {
    std::cout << "Filename was empty." << std::endl;
    return;
    }

  LoadPoints(this->LeftPane, fileName.toStdString());
}

void ResectioningWidget::on_actionLoadPointsRight_activated()
{
  // Get a filename to open
  QString fileName = QFileDialog::getOpenFileName(this, "Open File", ".", "Image Files (*.jpg *.jpeg *.bmp *.png *.mha)");

  std::cout << "Got filename: " << fileName.toStdString() << std::endl;
  if(fileName.toStdString().empty())
    {
    std::cout << "Filename was empty." << std::endl;
    return;
    }
  LoadPoints(this->RightPane, fileName.toStdString());
}

void ResectioningWidget::on_actionOpenImageLeft_activated()
{
  // Get a filename to open
  QString fileName = QFileDialog::getOpenFileName(this, "Open File", ".", "Image Files (*.jpg *.jpeg *.bmp *.png *.mha)");

  std::cout << "Got filename: " << fileName.toStdString() << std::endl;
  if(fileName.toStdString().empty())
    {
    std::cout << "Filename was empty." << std::endl;
    return;
    }
    
  if(this->LeftPane)
    {
    delete this->LeftPane;
    }
  this->LeftPane = new Pane2D(this->qvtkWidgetLeft);
  LoadImage(this->LeftPane, fileName.toStdString());
}

void ResectioningWidget::on_actionOpenImageRight_activated()
{
  // Get a filename to open
  QString fileName = QFileDialog::getOpenFileName(this, "Open File", ".", "Image Files (*.jpg *.jpeg *.bmp *.png *.mha)");

  std::cout << "Got filename: " << fileName.toStdString() << std::endl;
  if(fileName.toStdString().empty())
    {
    std::cout << "Filename was empty." << std::endl;
    return;
    }
    
  if(this->RightPane)
    {
    delete this->RightPane;
    }
  this->RightPane = new Pane2D(this->qvtkWidgetRight);
  LoadImage(this->RightPane, fileName.toStdString());
}


void ResectioningWidget::on_actionOpenPointCloudLeft_activated()
{
  // Get a filename to open
  QString fileName = QFileDialog::getOpenFileName(this, "Open Point Cloud", ".", "Files (*.vtp)");

  std::cout << "Got filename: " << fileName.toStdString() << std::endl;
  if(fileName.toStdString().empty())
    {
    std::cout << "Filename was empty." << std::endl;
    return;
    }
    
  if(this->LeftPane)
    {
    delete this->LeftPane;
    }
  this->LeftPane = new Pane3D(this->qvtkWidgetLeft);
  LoadPointCloud(this->LeftPane, fileName.toStdString());
  std::cout << "Done loading point cloud." << std::endl;
}

void ResectioningWidget::on_actionOpenPointCloudRight_activated()
{
  // Get a filename to open
  QString fileName = QFileDialog::getOpenFileName(this, "Open Point Cloud", ".", "Files (*.vtp)");

  std::cout << "Got filename: " << fileName.toStdString() << std::endl;
  if(fileName.toStdString().empty())
    {
    std::cout << "Filename was empty." << std::endl;
    return;
    }
    
  if(this->RightPane)
    {
    delete this->RightPane;
    }
  this->RightPane = new Pane3D(this->qvtkWidgetRight);
  LoadPointCloud(this->RightPane, fileName.toStdString());
  std::cout << "Done loading point cloud." << std::endl;
}

void ResectioningWidget::on_actionSavePointsLeft_activated()
{
  SavePoints(this->LeftPane);
}

void ResectioningWidget::on_actionSavePointsRight_activated()
{
  SavePoints(this->RightPane);
}

void ResectioningWidget::on_btnDeleteLastCorrespondenceLeft_clicked()
{
  this->LeftPane->SelectionStyle->DeleteLastCorrespondence();
  this->LeftPane->Refresh();
}

void ResectioningWidget::on_btnDeleteAllCorrespondencesLeft_clicked()
{
  this->LeftPane->SelectionStyle->RemoveAll();
  this->LeftPane->Refresh();
  //static_cast<PointSelectionStyle2D*>(pane->SelectionStyle)->RemoveAll();
  //this->qvtkWidgetLeft->GetRenderWindow()->Render();
}

void ResectioningWidget::on_btnDeleteLastCorrespondenceRight_clicked()
{
  this->RightPane->SelectionStyle->DeleteLastCorrespondence();
  this->RightPane->Refresh();
}

void ResectioningWidget::on_btnDeleteAllCorrespondencesRight_clicked()
{
  this->RightPane->SelectionStyle->RemoveAll();
  this->RightPane->Refresh();
  //static_cast<PointSelectionStyle2D*>(pane->SelectionStyle)->RemoveAll();
  //this->qvtkWidgetRight->GetRenderWindow()->Render();
}
