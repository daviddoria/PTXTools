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

#ifndef ResectioningWidget_H
#define ResectioningWidget_H

#include "ui_ResectioningWidget.h"

// ITK
#include "itkImage.h"

// Eigen
#include <Eigen/Dense>

// Qt
#include <QMainWindow>
#include <QFutureWatcher>
class QProgressDialog;

// VTK
class vtkSphereSource;

// Custom
#include "Pane.h"
#include "Pane2D.h"
#include "Pane3D.h"
#include "PointSelectionStyle.h"
#include "PTXImage.h"
#include "Types.h"

// Forward declarations
class vtkActor;
class vtkBorderWidget;
class vtkImageData;
class vtkImageActor;
class vtkPolyData;
class vtkPolyDataMapper;
class vtkRenderer;

class ResectioningWidget : public QMainWindow, public Ui::ResectioningWidget
{
  Q_OBJECT
public:

  /** Default constructor */
  ResectioningWidget();

  /** Constructor for loading files automatically. */
  ResectioningWidget(const std::string& imageFileName, const std::string& pointCloudFileName);

  /** Constructor for loading files and correspondences automatically. */
  ResectioningWidget(const std::string& imageFileName, const std::string& imageCorrespondenceFile,
                     const std::string& pointCloudFileName,
                     const std::string& pointCloudCorrespondenceFile);
  
  ~ResectioningWidget() {};

public slots:
  // Image
  void on_action_Image_Open_activated();
  void on_action_Image_SaveCorrespondences_activated();
  void on_action_Image_LoadCorrespondences_activated();
  void on_btnDeleteLastCorrespondenceImage_clicked();
  void on_btnDeleteAllCorrespondencesImage_clicked();
  void on_action_Image_FlipHorizontally_activated();
  void on_action_Image_FlipVertically_activated();

  // Point cloud
  void on_action_PointCloud_OpenVTP_activated();
  void on_action_PointCloud_OpenPTX_activated();
  void on_action_PointCloud_LoadCorrespondences_activated();
  void on_action_PointCloud_SaveCorrespondences_activated();
  void on_btnDeleteLastCorrespondencePointCloud_clicked();
  void on_btnDeleteAllCorrespondencesPointCloud_clicked();

  void on_chkInputMesh_clicked();
  void on_chkOutputMesh_clicked();
  
  // File menu
  void on_actionHelp_activated();
  void on_actionQuit_activated();

  // Export menu
  void on_action_Export_ResultPTX_activated();
  void on_action_Export_ResultRGB_activated();

  // Action buttons
  void on_btnResectionMesh_clicked();
  void on_btnResection_clicked();
  void on_btnComputeP_clicked();
  void on_btnWriteMesh_clicked();

private:

  Eigen::MatrixXd ComputeP();
    
  void SharedConstructor();

  /** Allow things to be run in a different thread while displaying a progress bar */
  QFutureWatcher<void> FutureWatcher;
  QProgressDialog* ProgressDialog;

  /** Load correspondences */
  void LoadCorrespondencesImage(const std::string& filename);
  void LoadCorrespondencesPointCloud(const std::string& filename);

  /** Save correspondences */
  void SaveCorrespondencesImage(const std::string& filename);
  void SaveCorrespondencesPointCloud(const std::string& filename);

  /** Load an image*/
  void LoadImage(const std::string& fileName);

  /** Get the resulting image from the resulting ptx and display it. */
  void ShowResultImage();

  /** Get the colored points from the resulting ptx and display them. */
  void ShowResultPointCloud();

  /** Load a VTP. */
  void LoadVTP(const std::string& fileName);

  /** Load a PTX. */
  void LoadPTX(const std::string& fileName);

  /** The pane responsible for displaying, manipulating, and selecting the image. */
  Pane2D* ImagePane;

  /** The pane responsible for displaying, manipulating, and selecting the point cloud. */
  Pane3D* PointCloudPane;

  /** The pane responsible for displaying, manipulating, and selecting the point cloud. */
  Pane3D* OutputPointCloudPane;
  
  /** The pane responsible for displaying the output image. */
  Pane2D* ResultImagePane;

  /** The active/displayed point cloud is produced from this. */
  PTXImage PTX;

  /** The point cloud that was loaded. We must save this because if a second resectioning is performed,
   *  it should start from the original, not the last iteration of resectioning. */
  PTXImage OriginalPTX;

  /** The image that was loaded. */
  FloatVectorImageType::Pointer ColorImage;
  //PTXImage::RGBImageType::Pointer ColorImage;

  vtkSmartPointer<vtkSphereSource> CameraLocationMarkerSource;
  vtkSmartPointer<vtkPolyDataMapper> CameraLocationMarkerMapper;
  vtkSmartPointer<vtkActor> CameraLocationMarker;

  vtkSmartPointer<vtkSphereSource> ScannerLocationMarkerSource;
  vtkSmartPointer<vtkPolyDataMapper> ScannerLocationMarkerMapper;
  vtkSmartPointer<vtkActor> ScannerLocationMarker;

  vtkSmartPointer<vtkPolyData> Mesh;

  // The "input" and "output" meshes are identical, but we need different actors for each renderer to follow good VTK practice
  vtkSmartPointer<vtkPolyDataMapper> InputMeshMapper;
  vtkSmartPointer<vtkActor> InputMeshActor;

  vtkSmartPointer<vtkPolyDataMapper> OutputMeshMapper;
  vtkSmartPointer<vtkActor> OutputMeshActor;
};

#endif
