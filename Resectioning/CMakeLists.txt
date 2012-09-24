cmake_minimum_required(VERSION 2.6)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

PROJECT(Resectioning)

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

if(${VTK_MAJOR_VERSION} LESS 6)
  message(FATAL_ERROR "Your VTK version (${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION}) is too old! At least 6.0 is required.")
endif()

# if(NOT VTK_USE_QVTK)
#   message(FATAL_ERROR "You must compile VTK with VTK_USE_QT=ON. To do this, you must also use VTK_USE_GUISUPPORT=ON and BUILD_SHARED_LIBS=ON.")
# endif()

if(NOT vtkGUISupportQt_LOADED)
  message(FATAL_ERROR "You must turn on the vtkGUISupportQt module in your VTK build.") # because the UI is Qt
endif()

if(NOT vtkFiltersFlowPaths_LOADED)
  message(FATAL_ERROR "You must turn on the vtkFiltersFlowPaths module in your VTK build.") # for vtkModifiedBSPTree
endif()

FIND_PACKAGE(Qt4 REQUIRED)
INCLUDE(${QT_USE_FILE})

FIND_PACKAGE(ITK 4 REQUIRED)
INCLUDE(${ITK_USE_FILE})

# We do this in the CameraCalibration submodule, but it doesn't seem to work unless it is here as well ("cannot find -lvgl" type of errors)
FIND_PACKAGE(VXL REQUIRED)
INCLUDE(${VXL_CMAKE_DIR}/UseVXL.cmake)

QT4_WRAP_UI(ResectioningUISrcs ResectioningWidget.ui)
QT4_WRAP_CPP(ResectioningMOCSrcs ResectioningWidget.h)

add_subdirectory(CameraCalibration)

ADD_EXECUTABLE(ResectioningInteractive
ResectioningInteractive.cpp
ResectioningWidget.cpp

Resectioning.cpp

ResectioningHelpers.cpp

Pane.cpp
Pane2D.cpp
Pane3D.cpp

Custom3DStyle.cpp

PointSelectionStyle.cpp
PointSelectionStyle2D.cpp
PointSelectionStyle3D.cpp

${ResectioningUISrcs} ${ResectioningMOCSrcs} ${ResourceSrcs} )
TARGET_LINK_LIBRARIES(ResectioningInteractive libCameraCalibration ${PTXTools_libraries})


# SET(ResectionSmart ON CACHE BOOL "Append one ptx to the right of another")
# if(ResectionSmart)
#   ADD_EXECUTABLE(ResectionSmart ColorPTXFromImage.cpp Resectioning.cpp ResectioningHelpers.cpp)
#   TARGET_LINK_LIBRARIES(ResectionSmart ${PTXLibraries})
#   INSTALL( TARGETS ResectionSmart RUNTIME DESTINATION ${INSTALL_DIR} )
# endif()
