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

// VTK
#include <vtkSmartPointer.h>
#include <vtkSeedWidget.h>
#include <vtkPointHandleRepresentation2D.h>

// ITK
#include "itkImage.h"

// Qt
#include <QMainWindow>
#include <QFutureWatcher>
class QProgressDialog;

// Custom
#include "Pane.h"
#include "Pane2D.h"
#include "Pane3D.h"
#include "PointSelectionStyle.h"
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
  
  ~ResectioningWidget() {};

public slots:
  // Image
  void on_action_Image_Open_activated();
  void on_action_Image_SaveCorrespondences_activated();
  void on_action_Image_LoadCorrespondences_activated();

  void on_action_PointCloud_OpenVTP_activated();
  void on_action_PointCloud_OpenPTX_activated();
  void on_action_PointCloud_LoadCorrespondences_activated();
  void on_action_PointCloud_SaveCorrespondences_activated();
  
  void on_actionHelp_activated();
  void on_actionQuit_activated();
  
  void on_btnDeleteLastCorrespondencePointCloud_clicked();
  void on_btnDeleteAllCorrespondencesPointCloud_clicked();
  void on_btnDeleteLastCorrespondenceImage_clicked();
  void on_btnDeleteAllCorrespondencesImage_clicked();

  void on_action_Image_FlipHorizontally_activated();
  void on_action_Image_FlipVertically_activated();

private:

  void SharedConstructor();

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

  /** Load a VTP. */
  void LoadVTP(const std::string& fileName);

  /** Load a PTX. */
  void LoadPTX(const std::string& fileName);
  
  Pane2D* ImagePane;
  Pane3D* PointCloudPane;
  
};

#endif
