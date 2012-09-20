cmake_minimum_required(VERSION 2.6)

PROJECT(PTXTools)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=gnu++0x")

SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

######### Libraries #########
#ITK
FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})
set(PTXTools_libraries ${PTXTools_libraries} ${ITK_LIBRARIES})

#VTK
FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})
set(PTXTools_libraries ${PTXTools_libraries} ${VTK_LIBRARIES})

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


# Submodules
# Only include ITKVTKHelpers if it is not already available.
get_property(FoundITKVTKHelpers GLOBAL PROPERTY ITKVTKHelpersIncludeDirs SET)
if(NOT FoundITKVTKHelpers)
  add_subdirectory(ITKVTKHelpers)
endif()

get_property(ITKVTKHelpersIncludeDirs GLOBAL PROPERTY ITKVTKHelpersIncludeDirs)
set(PTXTools_include_dirs ${PTXTools_include_dirs} ${ITKVTKHelpersIncludeDirs})
get_property(ITKVTKHelpersLibs GLOBAL PROPERTY ITKVTKHelpersLibs)
set(PTXTools_libraries ${PTXTools_libraries} ${ITKVTKHelpersLibs})

# Only include VTKHelpers if it is not already available.
get_property(FoundVTKHelpers GLOBAL PROPERTY VTKHelpersIncludeDirs SET)
if(NOT FoundVTKHelpers)
  add_subdirectory(VTKHelpers)
endif()

get_property(VTKHelpersIncludeDirs GLOBAL PROPERTY VTKHelpersIncludeDirs)
set(PTXTools_include_dirs ${PTXTools_include_dirs} ${VTKHelpersIncludeDirs})
get_property(VTKHelpersLibs GLOBAL PROPERTY VTKHelpersLibs)
set(PTXTools_libraries ${PTXTools_libraries} ${VTKHelpersLibs})

# Give the compiler all required include directories
include_directories(${PTXTools_include_dirs})

# Create the main library that will be used in all of the exectuables
add_library(PTX PTXPixel.cpp PTXImage.cpp PTXReader.cpp)
target_link_libraries(PTX ${PTXTools_libraries})

# Add this library to the list of libraries to build all of the executables against
set(PTXTools_libraries ${PTXTools_libraries} PTX)

#############
SET(PTX_CoordinateImages ON CACHE BOOL "Write images of the coordinates")
set(TargetList ${TargetList} PTX_CoordinateImages)
if(PTX_CoordinateImages)
  ADD_EXECUTABLE(ptxCoordinateImages CoordinateImages.cpp)
  TARGET_LINK_LIBRARIES(ptxCoordinateImages ${PTXTools_libraries})
  INSTALL( TARGETS ptxCoordinateImages RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

SET(PTX_AppendRight ON CACHE BOOL "Append one ptx to the right of another")
set(TargetList ${TargetList} PTX_AppendRight)
if(PTX_AppendRight)
  ADD_EXECUTABLE(ptxAppendRight AppendRight.cpp)
  TARGET_LINK_LIBRARIES(ptxAppendRight ${PTXTools_libraries})
  #TARGET_LINK_LIBRARIES(ptxAppendRight PTX ${ITK_LIBRARIES} ${VTK_LIBRARIES} ${SubmoduleLibraries})
  INSTALL( TARGETS ptxAppendRight RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

SET(PTX_CropPTX ON CACHE BOOL "Crop a ptx file")
set(TargetList ${TargetList} PTX_CropPTX)
if(PTX_CropPTX)
  ADD_EXECUTABLE(ptxCropPTX CropPTX.cpp)
  TARGET_LINK_LIBRARIES(ptxCropPTX ${PTXTools_libraries})
  INSTALL( TARGETS ptxCropPTX RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

SET(PTX_DerivativeImages OFF CACHE BOOL "Write images of the derivatives")
set(TargetList ${TargetList} PTX_DerivativeImages)
if(PTX_DerivativeImages)
  ADD_EXECUTABLE(ptxDerivativeImages DerivativeImages.cpp)
  TARGET_LINK_LIBRARIES(ptxDerivativeImages ${PTXTools_libraries})
  INSTALL( TARGETS ptxDerivativeImages RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

SET(PTX_LaplacianImages ON CACHE BOOL "Write images of the Laplacian")
set(TargetList ${TargetList} PTX_LaplacianImages)
if(PTX_LaplacianImages)
  ADD_EXECUTABLE(ptxLaplacianImages LaplacianImages.cpp)
  TARGET_LINK_LIBRARIES(ptxLaplacianImages ${PTXTools_libraries})
  INSTALL( TARGETS ptxLaplacianImages RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

SET(PTX_ExtractHoleBoundaryPoints ON CACHE BOOL "Determine pixels on the boundaries of holes")
set(TargetList ${TargetList} PTX_ExtractHoleBoundaryPoints)
if(PTX_ExtractHoleBoundaryPoints)
  ADD_EXECUTABLE(ptxExtractHoleBoundaryPoints ExtractHoleBoundaryPoints.cpp)
  TARGET_LINK_LIBRARIES(ptxExtractHoleBoundaryPoints ${PTXTools_libraries})
  INSTALL( TARGETS ptxExtractHoleBoundaryPoints RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

# PTX to Red, Green, Blue
SET(PTX_PTXtoRGB ON CACHE BOOL "Extract the colors (RGB) from a ptx file")
set(TargetList ${TargetList} PTX_PTXtoRGB)
if(PTX_PTXtoRGB)
  ADD_EXECUTABLE(ptxPTXtoRGB PTXtoRGB.cpp)
  TARGET_LINK_LIBRARIES(ptxPTXtoRGB ${PTXTools_libraries})
  INSTALL( TARGETS ptxPTXtoRGB RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

# PTX to Red, Green, Blue, Depth
SET(PTX_PTXtoRGBD ON CACHE BOOL "Extract the colors and depths from a ptx file")
set(TargetList ${TargetList} PTX_PTXtoRGBD)
if(PTX_PTXtoRGBD)
  ADD_EXECUTABLE(ptxPTXtoRGBD PTXtoRGBD.cpp)
  TARGET_LINK_LIBRARIES(ptxPTXtoRGBD ${PTXTools_libraries})
  INSTALL( TARGETS ptxPTXtoRGBD RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

# PTX to Red, Green, Blue, Depth, Valid
SET(PTX_PTXtoRGBDV ON CACHE BOOL "Extract the colors and depths from a ptx file")
set(TargetList ${TargetList} PTX_PTXtoRGBDV)
if(PTX_PTXtoRGBDV)
  ADD_EXECUTABLE(ptxPTXtoRGBDV PTXtoRGBDV.cpp)
  TARGET_LINK_LIBRARIES(ptxPTXtoRGBDV ${PTXTools_libraries})
  INSTALL( TARGETS ptxPTXtoRGBDV RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

SET(PTX_PTXtoPointCloud ON CACHE BOOL "Extract the point coordinates")
set(TargetList ${TargetList} PTX_PTXtoPointCloud)
if(PTX_PTXtoPointCloud)
  ADD_EXECUTABLE(ptxPTXtoPointCloud PTXtoPointCloud.cpp)
  TARGET_LINK_LIBRARIES(ptxPTXtoPointCloud ${PTXTools_libraries})
  INSTALL( TARGETS ptxPTXtoPointCloud RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

SET(PTX_PTXtoMesh ON CACHE BOOL "Extract the point coordinates")
set(TargetList ${TargetList} PTX_PTXtoMesh)
if(PTX_PTXtoMesh)
  ADD_EXECUTABLE(ptxPTXtoMesh PTXtoMesh.cpp)
  TARGET_LINK_LIBRARIES(ptxPTXtoMesh ${PTXTools_libraries})
  INSTALL( TARGETS ptxPTXtoMesh RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

SET(PTX_PTXtoDepthImage ON CACHE BOOL "Write the ptx depth image")
set(TargetList ${TargetList} PTX_PTXtoDepthImage)
if(PTX_PTXtoDepthImage)
  ADD_EXECUTABLE(ptxPTXtoDepthImage PTXtoDepthImage.cpp)
  TARGET_LINK_LIBRARIES(ptxPTXtoDepthImage ${PTXTools_libraries})
  INSTALL( TARGETS ptxPTXtoDepthImage RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

SET(PTX_PTXtoIntensityImage ON CACHE BOOL "Write the ptx intensity image")
set(TargetList ${TargetList} PTX_PTXtoIntensityImage)
if(PTX_PTXtoIntensityImage)
  ADD_EXECUTABLE(ptxPTXtoIntensityImage PTXtoIntensityImage.cpp)
  TARGET_LINK_LIBRARIES(ptxPTXtoIntensityImage ${PTXTools_libraries})
  INSTALL( TARGETS ptxPTXtoIntensityImage RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

SET(PTX_PTXtoDepthThresholdMask ON CACHE BOOL "Create a mask from a depth threshold")
set(TargetList ${TargetList} PTX_PTXtoDepthThresholdMask)
if(PTX_PTXtoDepthThresholdMask)
  ADD_EXECUTABLE(ptxPTXtoDepthThresholdMask PTXtoDepthThresholdMask.cpp)
  TARGET_LINK_LIBRARIES(ptxPTXtoDepthThresholdMask ${PTXTools_libraries})
  INSTALL( TARGETS ptxPTXtoDepthThresholdMask RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

SET(PTX_ReplaceDepthImage ON CACHE BOOL "Replace the depths of points in a ptx file")
set(TargetList ${TargetList} PTX_ReplaceDepthImage)
if(PTX_ReplaceDepthImage)
  ADD_EXECUTABLE(ptxReplaceDepthImage PTXReplaceDepthImage.cpp)
  TARGET_LINK_LIBRARIES(ptxReplaceDepthImage ${PTXTools_libraries})
  INSTALL( TARGETS ptxReplaceDepthImage RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

SET(PTX_ReplaceValidityImage ON CACHE BOOL "Replace the depths of points in a ptx file")
set(TargetList ${TargetList} PTX_ReplaceValidityImage)
if(PTX_ReplaceValidityImage)
  ADD_EXECUTABLE(ptxReplaceValidityImage PTXReplaceValidityImage.cpp)
  TARGET_LINK_LIBRARIES(ptxReplaceValidityImage ${PTXTools_libraries})
  INSTALL( TARGETS ptxReplaceValidityImage RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

SET(PTX_ReplaceRGBD ON CACHE BOOL "Replace the depths and colors of points in a ptx file")
set(TargetList ${TargetList} PTX_ReplaceRGBD)
if(PTX_ReplaceRGBD)
  ADD_EXECUTABLE(ptxReplaceRGBD PTXReplaceRGBD.cpp)
  TARGET_LINK_LIBRARIES(ptxReplaceRGBD ${PTXTools_libraries})
  INSTALL( TARGETS ptxReplaceRGBD RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

SET(PTX_ReplaceRGB ON CACHE BOOL "Replace the colors of points in a ptx file")
set(TargetList ${TargetList} PTX_ReplaceRGB)
if(PTX_ReplaceRGB)
  ADD_EXECUTABLE(ptxReplaceRGB PTXReplaceRGB.cpp)
  TARGET_LINK_LIBRARIES(ptxReplaceRGB ${PTXTools_libraries})
  INSTALL( TARGETS ptxReplaceRGB RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

SET(PTX_ReplaceXYZ ON CACHE BOOL "Replace the coordinates of points in a ptx file")
set(TargetList ${TargetList} PTX_ReplaceXYZ)
if(PTX_ReplaceXYZ)
  ADD_EXECUTABLE(ptxReplaceXYZ PTXReplaceXYZ.cpp)
  TARGET_LINK_LIBRARIES(ptxReplaceXYZ ${PTXTools_libraries})
  INSTALL( TARGETS ptxReplaceXYZ RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

SET(PTX_CoordinateImages ON CACHE BOOL "Write images of the coordinates")
set(TargetList ${TargetList} PTX_CoordinateImages)
if(PTX_CoordinateImages)
  ADD_EXECUTABLE(ptxPTXtoEverything PTXtoEverything.cpp)
  TARGET_LINK_LIBRARIES(ptxPTXtoEverything ${PTXTools_libraries})
  INSTALL( TARGETS ptxPTXtoEverything RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

SET(PTX_PTXtoValidityMask ON CACHE BOOL "Write the validity mask")
set(TargetList ${TargetList} PTX_PTXtoValidityMask)
if(PTX_PTXtoValidityMask)
  ADD_EXECUTABLE(ptxPTXtoValidityMask PTXtoValidityMask.cpp)
  TARGET_LINK_LIBRARIES(ptxPTXtoValidityMask ${PTXTools_libraries})
  INSTALL( TARGETS ptxPTXtoValidityMask RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

SET(PTX_WeightedLaplacian ON CACHE BOOL "Compute and write the weight Laplacian")
set(TargetList ${TargetList} PTX_WeightedLaplacian)
if(PTX_WeightedLaplacian)
  ADD_EXECUTABLE(ptxWeightedLaplacian PTXWeightedLaplacian.cpp)
  TARGET_LINK_LIBRARIES(ptxWeightedLaplacian ${PTXTools_libraries})
  INSTALL( TARGETS ptxWeightedLaplacian RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

SET(PTX_Downsample ON CACHE BOOL "Downsample a ptx file")
set(TargetList ${TargetList} PTX_Downsample)
if(PTX_Downsample)
  ADD_EXECUTABLE(ptxDownsample Downsample.cpp)
  TARGET_LINK_LIBRARIES(ptxDownsample ${PTXTools_libraries})
  INSTALL( TARGETS ptxDownsample RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

# SET(PTX_ColorPTXFromImage ON CACHE BOOL "Remap the colors from an image to the ptx file")
# set(TargetList ${TargetList} PTX_ColorPTXFromImage)
# if(PTX_ColorPTXFromImage)
#   ADD_EXECUTABLE(ptxColorPTXFromImage ColorPTXFromImage.cpp)
#   TARGET_LINK_LIBRARIES(ptxColorPTXFromImage ${PTXTools_libraries})
#   INSTALL( TARGETS ptxColorPTXFromImage RUNTIME DESTINATION ${INSTALL_DIR} )
# endif()
#
# SET(PTX_ColorPTXFromImageNaive ON CACHE BOOL "Remap the colors from an image to the ptx file")
# set(TargetList ${TargetList} PTX_ColorPTXFromImageNaive)
# if(PTX_ColorPTXFromImageNaive)
#   ADD_EXECUTABLE(ptxColorPTXFromImageNaive ColorPTXFromImageNaive.cpp)
#   TARGET_LINK_LIBRARIES(ptxColorPTXFromImageNaive ${PTXTools_libraries})
#   INSTALL( TARGETS ptxColorPTXFromImageNaive RUNTIME DESTINATION ${INSTALL_DIR} )
# endif()

SET(PTX_ExtractMaskedRegion ON CACHE BOOL "Extract a part of a ptx file")
set(TargetList ${TargetList} PTX_ExtractMaskedRegion)
if(PTX_ExtractMaskedRegion)
  ADD_EXECUTABLE(ptxExtractMaskedRegion ExtractMaskedRegion.cpp)
  TARGET_LINK_LIBRARIES(ptxExtractMaskedRegion ${PTXTools_libraries})
  INSTALL( TARGETS ptxExtractMaskedRegion RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

SET(PTX_OrthogonalProjection ON CACHE BOOL "Orthogonally project the points in a ptx file")
set(TargetList ${TargetList} PTX_OrthogonalProjection)
if(PTX_OrthogonalProjection)
  ADD_EXECUTABLE(ptxOrthogonalProjection OrthogonalProjection.cpp)
  TARGET_LINK_LIBRARIES(ptxOrthogonalProjection ${PTXTools_libraries})
  set_target_properties(ptxOrthogonalProjection PROPERTIES COMPILE_DEFINITIONS "${DEFINITIONS}")
  INSTALL( TARGETS ptxOrthogonalProjection RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

SET(PTX_PTXInfo ON CACHE BOOL "Orthogonally project the points in a ptx file")
set(TargetList ${TargetList} PTX_PTXInfo)
if(PTX_PTXInfo)
  ADD_EXECUTABLE(ptxPTXInfo PTXInfo.cpp)
  TARGET_LINK_LIBRARIES(ptxPTXInfo ${PTXTools_libraries})
  set_target_properties(ptxPTXInfo PROPERTIES COMPILE_DEFINITIONS "${DEFINITIONS}")
  INSTALL( TARGETS ptxPTXInfo RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

SET(PTX_PointCloudPropertyToImage ON CACHE BOOL "Create an image from a property attached to points in a point cloud")
set(TargetList ${TargetList} PTX_PointCloudPropertyToImage)
if(PTX_PointCloudPropertyToImage)
  ADD_EXECUTABLE(ptxPointCloudPropertyToImage PointCloudPropertyToImage.cpp)
  TARGET_LINK_LIBRARIES(ptxPointCloudPropertyToImage ${PTXTools_libraries})
  set_target_properties(ptxPointCloudPropertyToImage PROPERTIES COMPILE_DEFINITIONS "${DEFINITIONS}")
  INSTALL( TARGETS ptxPointCloudPropertyToImage RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

SET(PTX_MakeAllPointsValid ON CACHE BOOL "Mark all points in the PTX as valid.")
set(TargetList ${TargetList} PTX_MakeAllPointsValid)
if(PTX_MakeAllPointsValid)
  ADD_EXECUTABLE(ptxMakeAllPointsValid PTXMakeAllPointsValid.cpp)
  TARGET_LINK_LIBRARIES(ptxMakeAllPointsValid ${PTXTools_libraries})
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

# Allow this project to be detected and used as a submodule
set(PTXTools_include_dirs ${PTXTools_include_dirs} ${CMAKE_CURRENT_SOURCE_DIR})
set_property(GLOBAL PROPERTY PTXToolsIncludeDirs ${PTXTools_include_dirs})
set_property(GLOBAL PROPERTY PTXToolsLibs ${PTXTools_libraries})

####### Mass toggle options #######
option(TurnOnAllTargets "all")
option(TurnOffAllTargets "all")

if(TurnOnAllTargets)
  foreach(CurrentTarget ${TargetList})
    set(${CurrentTarget} ON CACHE BOOL "test" FORCE)
  endforeach()
  set(TurnOnAllTargets OFF CACHE BOOL "test" FORCE)  # Set itself back to off, as this is a one time thing.
endif()

if(TurnOffAllTargets)
  foreach(CurrentTarget ${TargetList})
    set(${CurrentTarget} OFF CACHE BOOL "test" FORCE)
  endforeach()
  set(TurnOffAllTargets OFF CACHE BOOL "test" FORCE) # Set itself back to off, as this is a one time thing.
endif()
