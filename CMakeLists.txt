cmake_minimum_required(VERSION 2.6)

PROJECT(PTXTools)
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

######### Libraries #########3
#ITK
FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})

#VTK
FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

#Eigen
#set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR} ${CMAKE_MODULE_PATH})
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR})

# FIND_PACKAGE(Eigen3 REQUIRED)
# INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})

# Need to set EIGEN3_INCLUDE_DIR in this project CMake settings (or do the below:)
include_directories(/media/portable/src/Eigen)

# Setup install
SET( INSTALL_DIR ${CMAKE_INSTALL_PREFIX} )

# if(CMAKE_BUILD_TYPE STREQUAL "Release")
#   message("Compiled in Release mode")
# else(CMAKE_BUILD_TYPE STREQUAL "Release")
#   list(APPEND DEFINITIONS "DEBUGMODE")
#   message("Compiled in Debug mode")
# endif(CMAKE_BUILD_TYPE STREQUAL "Release")

# Create the main library that will be used in all of the exectuables
add_library(PTX PTXPixel.cpp PTXImage.cpp PTXReader.cpp Helpers.cpp)
set(PTXLibraries PTX ${ITK_LIBRARIES} ${VTK_LIBRARIES})

#############
SET(PTX_CoordinateImages ON CACHE BOOL "Write images of the coordinates")
if(PTX_CoordinateImages)
  ADD_EXECUTABLE(ptxCoordinateImages CoordinateImages.cpp)
  TARGET_LINK_LIBRARIES(ptxCoordinateImages ${PTXLibraries})
  INSTALL( TARGETS ptxCoordinateImages RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

SET(PTX_AppendRight ON CACHE BOOL "Append one ptx to the right of another")
if(PTX_AppendRight)
  ADD_EXECUTABLE(ptxAppendRight AppendRight.cpp)
  TARGET_LINK_LIBRARIES(ptxAppendRight ${PTXLibraries})
  INSTALL( TARGETS ptxAppendRight RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

SET(PTX_CropPTX ON CACHE BOOL "Crop a ptx file")
if(PTX_CropPTX)
  ADD_EXECUTABLE(ptxCropPTX CropPTX.cpp)
  TARGET_LINK_LIBRARIES(ptxCropPTX ${PTXLibraries})
  INSTALL( TARGETS ptxCropPTX RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

SET(PTX_DerivativeImages OFF CACHE BOOL "Write images of the derivatives")
if(PTX_DerivativeImages)
  ADD_EXECUTABLE(ptxDerivativeImages DerivativeImages.cpp)
  TARGET_LINK_LIBRARIES(ptxDerivativeImages ${PTXLibraries})
  INSTALL( TARGETS ptxDerivativeImages RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

SET(PTX_LaplacianImages ON CACHE BOOL "Write images of the Laplacian")
if(PTX_LaplacianImages)
  ADD_EXECUTABLE(ptxLaplacianImages LaplacianImages.cpp)
  TARGET_LINK_LIBRARIES(ptxLaplacianImages ${PTXLibraries})
  INSTALL( TARGETS ptxLaplacianImages RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

SET(PTX_ExtractHoleBoundaryPoints ON CACHE BOOL "Determine pixels on the boundaries of holes")
if(PTX_ExtractHoleBoundaryPoints)
  ADD_EXECUTABLE(ptxExtractHoleBoundaryPoints ExtractHoleBoundaryPoints.cpp)
  TARGET_LINK_LIBRARIES(ptxExtractHoleBoundaryPoints ${PTXLibraries})
  INSTALL( TARGETS ptxExtractHoleBoundaryPoints RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

# PTX to Red, Green, Blue
SET(PTX_PTXtoRGB ON CACHE BOOL "Extract the colors (RGB) from a ptx file")
if(PTX_PTXtoRGB)
  ADD_EXECUTABLE(ptxPTXtoRGB PTXtoRGB.cpp)
  TARGET_LINK_LIBRARIES(ptxPTXtoRGB ${PTXLibraries})
  INSTALL( TARGETS ptxPTXtoRGB RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

# PTX to Red, Green, Blue, Depth
SET(PTX_PTXtoRGBD ON CACHE BOOL "Extract the colors and depths from a ptx file")
if(PTX_PTXtoRGBD)
  ADD_EXECUTABLE(ptxPTXtoRGBD PTXtoRGBD.cpp)
  TARGET_LINK_LIBRARIES(ptxPTXtoRGBD ${PTXLibraries})
  INSTALL( TARGETS ptxPTXtoRGBD RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

# PTX to Red, Green, Blue, Depth, Valid
SET(PTX_PTXtoRGBDV ON CACHE BOOL "Extract the colors and depths from a ptx file")
if(PTX_PTXtoRGBDV)
  ADD_EXECUTABLE(ptxPTXtoRGBDV PTXtoRGBDV.cpp)
  TARGET_LINK_LIBRARIES(ptxPTXtoRGBDV ${PTXLibraries})
  INSTALL( TARGETS ptxPTXtoRGBDV RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

SET(PTX_PTXtoPointCloud ON CACHE BOOL "Extract the point coordinates")
if(PTX_PTXtoPointCloud)
  ADD_EXECUTABLE(ptxPTXtoPointCloud PTXtoPointCloud.cpp)
  TARGET_LINK_LIBRARIES(ptxPTXtoPointCloud ${PTXLibraries})
  INSTALL( TARGETS ptxPTXtoPointCloud RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

SET(PTX_PTXtoMesh ON CACHE BOOL "Extract the point coordinates")
if(PTX_PTXtoMesh)
  ADD_EXECUTABLE(ptxPTXtoMesh PTXtoMesh.cpp)
  TARGET_LINK_LIBRARIES(ptxPTXtoMesh ${PTXLibraries})
  INSTALL( TARGETS ptxPTXtoMesh RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

SET(PTX_PTXtoDepthImage ON CACHE BOOL "Write the ptx depth image")
if(PTX_PTXtoDepthImage)
  ADD_EXECUTABLE(ptxPTXtoDepthImage PTXtoDepthImage.cpp)
  TARGET_LINK_LIBRARIES(ptxPTXtoDepthImage ${PTXLibraries})
  INSTALL( TARGETS ptxPTXtoDepthImage RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

SET(PTX_PTXtoIntensityImage ON CACHE BOOL "Write the ptx intensity image")
if(PTX_PTXtoIntensityImage)
  ADD_EXECUTABLE(ptxPTXtoIntensityImage PTXtoIntensityImage.cpp)
  TARGET_LINK_LIBRARIES(ptxPTXtoIntensityImage ${PTXLibraries})
  INSTALL( TARGETS ptxPTXtoIntensityImage RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

SET(PTX_PTXtoDepthThresholdMask ON CACHE BOOL "Create a mask from a depth threshold")
if(PTX_PTXtoDepthThresholdMask)
  ADD_EXECUTABLE(ptxPTXtoDepthThresholdMask PTXtoDepthThresholdMask.cpp)
  TARGET_LINK_LIBRARIES(ptxPTXtoDepthThresholdMask ${PTXLibraries})
  INSTALL( TARGETS ptxPTXtoDepthThresholdMask RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

SET(PTX_ReplaceDepthImage ON CACHE BOOL "Replace the depths of points in a ptx file")
if(PTX_ReplaceDepthImage)
  ADD_EXECUTABLE(ptxReplaceDepthImage PTXReplaceDepthImage.cpp)
  TARGET_LINK_LIBRARIES(ptxReplaceDepthImage ${PTXLibraries})
  INSTALL( TARGETS ptxReplaceDepthImage RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

SET(PTX_ReplaceValidityImage ON CACHE BOOL "Replace the depths of points in a ptx file")
if(PTX_ReplaceValidityImage)
  ADD_EXECUTABLE(ptxReplaceValidityImage PTXReplaceValidityImage.cpp)
  TARGET_LINK_LIBRARIES(ptxReplaceValidityImage ${PTXLibraries})
  INSTALL( TARGETS ptxReplaceValidityImage RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

SET(PTX_ReplaceRGBD ON CACHE BOOL "Replace the depths and colors of points in a ptx file")
if(PTX_ReplaceRGBD)
  ADD_EXECUTABLE(ptxReplaceRGBD PTXReplaceRGBD.cpp)
  TARGET_LINK_LIBRARIES(ptxReplaceRGBD ${PTXLibraries})
  INSTALL( TARGETS ptxReplaceRGBD RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

SET(PTX_ReplaceRGB ON CACHE BOOL "Replace the colors of points in a ptx file")
if(PTX_ReplaceRGB)
  ADD_EXECUTABLE(ptxReplaceRGB PTXReplaceRGB.cpp)
  TARGET_LINK_LIBRARIES(ptxReplaceRGB ${PTXLibraries})
  INSTALL( TARGETS ptxReplaceRGB RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

SET(PTX_ReplaceXYZ ON CACHE BOOL "Replace the coordinates of points in a ptx file")
if(PTX_ReplaceXYZ)
  ADD_EXECUTABLE(ptxReplaceXYZ PTXReplaceXYZ.cpp)
  TARGET_LINK_LIBRARIES(ptxReplaceXYZ ${PTXLibraries})
  INSTALL( TARGETS ptxReplaceXYZ RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

SET(PTX_CoordinateImages ON CACHE BOOL "Write images of the coordinates")
if(PTX_CoordinateImages)
  ADD_EXECUTABLE(ptxPTXtoEverything PTXtoEverything.cpp)
  TARGET_LINK_LIBRARIES(ptxPTXtoEverything ${PTXLibraries})
  INSTALL( TARGETS ptxPTXtoEverything RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

SET(PTX_PTXtoValidityMask ON CACHE BOOL "Write the validity mask")
if(PTX_PTXtoValidityMask)
  ADD_EXECUTABLE(ptxPTXtoValidityMask PTXtoValidityMask.cpp)
  TARGET_LINK_LIBRARIES(ptxPTXtoValidityMask ${PTXLibraries})
  INSTALL( TARGETS ptxPTXtoValidityMask RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

SET(PTX_WeightedLaplacian ON CACHE BOOL "Compute and write the weight Laplacian")
if(PTX_WeightedLaplacian)
  ADD_EXECUTABLE(ptxWeightedLaplacian PTXWeightedLaplacian.cpp)
  TARGET_LINK_LIBRARIES(ptxWeightedLaplacian ${PTXLibraries})
  INSTALL( TARGETS ptxWeightedLaplacian RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

SET(PTX_Downsample ON CACHE BOOL "Downsample a ptx file")
if(PTX_Downsample)
  ADD_EXECUTABLE(ptxDownsample Downsample.cpp)
  TARGET_LINK_LIBRARIES(ptxDownsample ${PTXLibraries})
  INSTALL( TARGETS ptxDownsample RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

# SET(PTX_ColorPTXFromImage ON CACHE BOOL "Remap the colors from an image to the ptx file")
# if(PTX_ColorPTXFromImage)
#   ADD_EXECUTABLE(ptxColorPTXFromImage ColorPTXFromImage.cpp)
#   TARGET_LINK_LIBRARIES(ptxColorPTXFromImage ${PTXLibraries})
#   INSTALL( TARGETS ptxColorPTXFromImage RUNTIME DESTINATION ${INSTALL_DIR} )
# endif()
# 
# SET(PTX_ColorPTXFromImageNaive ON CACHE BOOL "Remap the colors from an image to the ptx file")
# if(PTX_ColorPTXFromImageNaive)
#   ADD_EXECUTABLE(ptxColorPTXFromImageNaive ColorPTXFromImageNaive.cpp)
#   TARGET_LINK_LIBRARIES(ptxColorPTXFromImageNaive ${PTXLibraries})
#   INSTALL( TARGETS ptxColorPTXFromImageNaive RUNTIME DESTINATION ${INSTALL_DIR} )
# endif()

SET(PTX_ExtractMaskedRegion ON CACHE BOOL "Extract a part of a ptx file")
if(PTX_ExtractMaskedRegion)
  ADD_EXECUTABLE(ptxExtractMaskedRegion ExtractMaskedRegion.cpp)
  TARGET_LINK_LIBRARIES(ptxExtractMaskedRegion ${PTXLibraries})
  INSTALL( TARGETS ptxExtractMaskedRegion RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

SET(PTX_OrthogonalProjection ON CACHE BOOL "Orthogonally project the points in a ptx file")
if(PTX_OrthogonalProjection)
  ADD_EXECUTABLE(ptxOrthogonalProjection OrthogonalProjection.cpp)
  TARGET_LINK_LIBRARIES(ptxOrthogonalProjection ${PTXLibraries})
  set_target_properties(ptxOrthogonalProjection PROPERTIES COMPILE_DEFINITIONS "${DEFINITIONS}")
  INSTALL( TARGETS ptxOrthogonalProjection RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

SET(PTX_PTXInfo ON CACHE BOOL "Orthogonally project the points in a ptx file")
if(PTX_PTXInfo)
  ADD_EXECUTABLE(ptxPTXInfo PTXInfo.cpp)
  TARGET_LINK_LIBRARIES(ptxPTXInfo ${PTXLibraries})
  set_target_properties(ptxPTXInfo PROPERTIES COMPILE_DEFINITIONS "${DEFINITIONS}")
  INSTALL( TARGETS ptxPTXInfo RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

SET(PTX_PointCloudPropertyToImage ON CACHE BOOL "Create an image from a property attached to points in a point cloud")
if(PTX_PointCloudPropertyToImage)
  ADD_EXECUTABLE(ptxPointCloudPropertyToImage PointCloudPropertyToImage.cpp)
  TARGET_LINK_LIBRARIES(ptxPointCloudPropertyToImage ${PTXLibraries})
  set_target_properties(ptxPointCloudPropertyToImage PROPERTIES COMPILE_DEFINITIONS "${DEFINITIONS}")
  INSTALL( TARGETS ptxPointCloudPropertyToImage RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

SET(PTX_MakeAllPointsValid ON CACHE BOOL "Mark all points in the PTX as valid.")
if(PTX_MakeAllPointsValid)
  ADD_EXECUTABLE(ptxMakeAllPointsValid PTXMakeAllPointsValid.cpp)
  TARGET_LINK_LIBRARIES(ptxMakeAllPointsValid ${PTXLibraries})
  set_target_properties(ptxMakeAllPointsValid PROPERTIES COMPILE_DEFINITIONS "${DEFINITIONS}")
  INSTALL( TARGETS ptxMakeAllPointsValid RUNTIME DESTINATION ${INSTALL_DIR} )
endif()


#############################################
SET(PTX_PTXViewer ON CACHE BOOL "An interactive ptx file viewer.")
if(PTX_PTXViewer)
  add_subdirectory(GUI)
endif()

SET(PTX_Resectioning ON CACHE BOOL "Color the PTX by an external image")
if(PTX_Resectioning)
  add_subdirectory(Resectioning)
endif()