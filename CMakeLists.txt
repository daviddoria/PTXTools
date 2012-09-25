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
option(PTX_CoordinateImages "Write images of the coordinates" ON)
if(PTX_CoordinateImages)
  ADD_EXECUTABLE(ptxCoordinateImages CoordinateImages.cpp)
  TARGET_LINK_LIBRARIES(ptxCoordinateImages ${PTXTools_libraries})
  INSTALL( TARGETS ptxCoordinateImages RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

option(PTX_AppendRight "Append one ptx to the right of another" ON)
if(PTX_AppendRight)
  ADD_EXECUTABLE(ptxAppendRight AppendRight.cpp)
  TARGET_LINK_LIBRARIES(ptxAppendRight ${PTXTools_libraries})
  #TARGET_LINK_LIBRARIES(ptxAppendRight PTX ${ITK_LIBRARIES} ${VTK_LIBRARIES} ${SubmoduleLibraries})
  INSTALL( TARGETS ptxAppendRight RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

option(PTX_CropPTX "Crop a ptx file" ON)
if(PTX_CropPTX)
  ADD_EXECUTABLE(ptxCropPTX CropPTX.cpp)
  TARGET_LINK_LIBRARIES(ptxCropPTX ${PTXTools_libraries})
  INSTALL( TARGETS ptxCropPTX RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

option(PTX_DerivativeImages "Write images of the derivatives" ON)
if(PTX_DerivativeImages)
  ADD_EXECUTABLE(ptxDerivativeImages DerivativeImages.cpp)
  TARGET_LINK_LIBRARIES(ptxDerivativeImages ${PTXTools_libraries})
  INSTALL( TARGETS ptxDerivativeImages RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

option(PTX_LaplacianImages "Write images of the Laplacian" ON)
if(PTX_LaplacianImages)
  ADD_EXECUTABLE(ptxLaplacianImages LaplacianImages.cpp)
  TARGET_LINK_LIBRARIES(ptxLaplacianImages ${PTXTools_libraries})
  INSTALL( TARGETS ptxLaplacianImages RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

option(PTX_ExtractHoleBoundaryPoints "Determine pixels on the boundaries of holes" ON)
if(PTX_ExtractHoleBoundaryPoints)
  ADD_EXECUTABLE(ptxExtractHoleBoundaryPoints ExtractHoleBoundaryPoints.cpp)
  TARGET_LINK_LIBRARIES(ptxExtractHoleBoundaryPoints ${PTXTools_libraries})
  INSTALL( TARGETS ptxExtractHoleBoundaryPoints RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

# PTX to Red, Green, Blue
option(PTX_PTXtoRGB "Extract the colors (RGB) from a ptx file" ON)
if(PTX_PTXtoRGB)
  ADD_EXECUTABLE(ptxPTXtoRGB PTXtoRGB.cpp)
  TARGET_LINK_LIBRARIES(ptxPTXtoRGB ${PTXTools_libraries})
  INSTALL( TARGETS ptxPTXtoRGB RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

# PTX to Red, Green, Blue, Depth
option(PTX_PTXtoRGBD "Extract the colors and depths from a ptx file" ON)
if(PTX_PTXtoRGBD)
  ADD_EXECUTABLE(ptxPTXtoRGBD PTXtoRGBD.cpp)
  TARGET_LINK_LIBRARIES(ptxPTXtoRGBD ${PTXTools_libraries})
  INSTALL( TARGETS ptxPTXtoRGBD RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

# PTX to Red, Green, Blue, Depth, Valid
option(PTX_PTXtoRGBDV "Extract the colors and depths from a ptx file" ON)
if(PTX_PTXtoRGBDV)
  ADD_EXECUTABLE(ptxPTXtoRGBDV PTXtoRGBDV.cpp)
  TARGET_LINK_LIBRARIES(ptxPTXtoRGBDV ${PTXTools_libraries})
  INSTALL( TARGETS ptxPTXtoRGBDV RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

option(PTX_PTXtoPointCloud "Extract the point coordinates" ON)
if(PTX_PTXtoPointCloud)
  ADD_EXECUTABLE(ptxPTXtoPointCloud PTXtoPointCloud.cpp)
  TARGET_LINK_LIBRARIES(ptxPTXtoPointCloud ${PTXTools_libraries})
  INSTALL( TARGETS ptxPTXtoPointCloud RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

option(PTX_PTXtoMesh "Extract the point coordinates" ON)
if(PTX_PTXtoMesh)
  ADD_EXECUTABLE(ptxPTXtoMesh PTXtoMesh.cpp)
  TARGET_LINK_LIBRARIES(ptxPTXtoMesh ${PTXTools_libraries})
  INSTALL( TARGETS ptxPTXtoMesh RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

option(PTX_PTXtoDepthImage "Write the ptx depth image" ON)
if(PTX_PTXtoDepthImage)
  ADD_EXECUTABLE(ptxPTXtoDepthImage PTXtoDepthImage.cpp)
  TARGET_LINK_LIBRARIES(ptxPTXtoDepthImage ${PTXTools_libraries})
  INSTALL( TARGETS ptxPTXtoDepthImage RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

option(PTX_PTXtoIntensityImage "Write the ptx intensity image" ON)
if(PTX_PTXtoIntensityImage)
  ADD_EXECUTABLE(ptxPTXtoIntensityImage PTXtoIntensityImage.cpp)
  TARGET_LINK_LIBRARIES(ptxPTXtoIntensityImage ${PTXTools_libraries})
  INSTALL( TARGETS ptxPTXtoIntensityImage RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

option(PTX_PTXtoDepthThresholdMask "Create a mask from a depth threshold" ON)
if(PTX_PTXtoDepthThresholdMask)
  ADD_EXECUTABLE(ptxPTXtoDepthThresholdMask PTXtoDepthThresholdMask.cpp)
  TARGET_LINK_LIBRARIES(ptxPTXtoDepthThresholdMask ${PTXTools_libraries})
  INSTALL( TARGETS ptxPTXtoDepthThresholdMask RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

option(PTX_ReplaceDepthImage "Replace the depths of points in a ptx file" ON)
if(PTX_ReplaceDepthImage)
  ADD_EXECUTABLE(ptxReplaceDepthImage PTXReplaceDepthImage.cpp)
  TARGET_LINK_LIBRARIES(ptxReplaceDepthImage ${PTXTools_libraries})
  INSTALL( TARGETS ptxReplaceDepthImage RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

option(PTX_ReplaceValidityImage "Replace the depths of points in a ptx file" ON)
if(PTX_ReplaceValidityImage)
  ADD_EXECUTABLE(ptxReplaceValidityImage PTXReplaceValidityImage.cpp)
  TARGET_LINK_LIBRARIES(ptxReplaceValidityImage ${PTXTools_libraries})
  INSTALL( TARGETS ptxReplaceValidityImage RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

option(PTX_ReplaceRGBD "Replace the depths and colors of points in a ptx file" ON)
if(PTX_ReplaceRGBD)
  ADD_EXECUTABLE(ptxReplaceRGBD PTXReplaceRGBD.cpp)
  TARGET_LINK_LIBRARIES(ptxReplaceRGBD ${PTXTools_libraries})
  INSTALL( TARGETS ptxReplaceRGBD RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

option(PTX_ReplaceRGB "Replace the colors of points in a ptx file" ON)
if(PTX_ReplaceRGB)
  ADD_EXECUTABLE(ptxReplaceRGB PTXReplaceRGB.cpp)
  TARGET_LINK_LIBRARIES(ptxReplaceRGB ${PTXTools_libraries})
  INSTALL( TARGETS ptxReplaceRGB RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

option(PTX_ReplaceXYZ "Replace the coordinates of points in a ptx file" ON)
if(PTX_ReplaceXYZ)
  ADD_EXECUTABLE(ptxReplaceXYZ PTXReplaceXYZ.cpp)
  TARGET_LINK_LIBRARIES(ptxReplaceXYZ ${PTXTools_libraries})
  INSTALL( TARGETS ptxReplaceXYZ RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

option(PTX_CoordinateImages "Write images of the coordinates" ON)
if(PTX_CoordinateImages)
  ADD_EXECUTABLE(ptxPTXtoEverything PTXtoEverything.cpp)
  TARGET_LINK_LIBRARIES(ptxPTXtoEverything ${PTXTools_libraries})
  INSTALL( TARGETS ptxPTXtoEverything RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

option(PTX_PTXtoValidityMask "Write the validity mask" ON)
if(PTX_PTXtoValidityMask)
  ADD_EXECUTABLE(ptxPTXtoValidityMask PTXtoValidityMask.cpp)
  TARGET_LINK_LIBRARIES(ptxPTXtoValidityMask ${PTXTools_libraries})
  INSTALL( TARGETS ptxPTXtoValidityMask RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

option(PTX_WeightedLaplacian "Compute and write the weight Laplacian" ON)
if(PTX_WeightedLaplacian)
  ADD_EXECUTABLE(ptxWeightedLaplacian PTXWeightedLaplacian.cpp)
  TARGET_LINK_LIBRARIES(ptxWeightedLaplacian ${PTXTools_libraries})
  INSTALL( TARGETS ptxWeightedLaplacian RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

option(PTX_Downsample "Downsample a ptx file" ON)
if(PTX_Downsample)
  ADD_EXECUTABLE(ptxDownsample Downsample.cpp)
  TARGET_LINK_LIBRARIES(ptxDownsample ${PTXTools_libraries})
  INSTALL( TARGETS ptxDownsample RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

option(PTX_ExtractMaskedRegion "Extract a part of a ptx file" ON)
if(PTX_ExtractMaskedRegion)
  ADD_EXECUTABLE(ptxExtractMaskedRegion ExtractMaskedRegion.cpp)
  TARGET_LINK_LIBRARIES(ptxExtractMaskedRegion ${PTXTools_libraries})
  INSTALL( TARGETS ptxExtractMaskedRegion RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

option(PTX_OrthogonalProjection "Orthogonally project the points in a ptx file" ON)
if(PTX_OrthogonalProjection)
  ADD_EXECUTABLE(ptxOrthogonalProjection OrthogonalProjection.cpp)
  TARGET_LINK_LIBRARIES(ptxOrthogonalProjection ${PTXTools_libraries})
  set_target_properties(ptxOrthogonalProjection PROPERTIES COMPILE_DEFINITIONS "${DEFINITIONS}")
  INSTALL( TARGETS ptxOrthogonalProjection RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

option(PTX_PTXInfo "Orthogonally project the points in a ptx file" ON)
if(PTX_PTXInfo)
  ADD_EXECUTABLE(ptxPTXInfo PTXInfo.cpp)
  TARGET_LINK_LIBRARIES(ptxPTXInfo ${PTXTools_libraries})
  set_target_properties(ptxPTXInfo PROPERTIES COMPILE_DEFINITIONS "${DEFINITIONS}")
  INSTALL( TARGETS ptxPTXInfo RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

option(PTX_PointCloudPropertyToImage "Create an image from a property attached to points in a point cloud" ON)
if(PTX_PointCloudPropertyToImage)
  ADD_EXECUTABLE(ptxPointCloudPropertyToImage PointCloudPropertyToImage.cpp)
  TARGET_LINK_LIBRARIES(ptxPointCloudPropertyToImage ${PTXTools_libraries})
  set_target_properties(ptxPointCloudPropertyToImage PROPERTIES COMPILE_DEFINITIONS "${DEFINITIONS}")
  INSTALL( TARGETS ptxPointCloudPropertyToImage RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

option(PTX_MakeAllPointsValid "Mark all points in the PTX as valid." ON)
if(PTX_MakeAllPointsValid)
  ADD_EXECUTABLE(ptxMakeAllPointsValid PTXMakeAllPointsValid.cpp)
  TARGET_LINK_LIBRARIES(ptxMakeAllPointsValid ${PTXTools_libraries})
  set_target_properties(ptxMakeAllPointsValid PROPERTIES COMPILE_DEFINITIONS "${DEFINITIONS}")
  INSTALL( TARGETS ptxMakeAllPointsValid RUNTIME DESTINATION ${INSTALL_DIR} )
endif()

#############################################
# These variables do not begin with PTX_ because they are not logically in the same "group" as the "small PTX executables" above,
# so we don't wanto to enable them unless the user explicitly requests them.
option(PTXViewer "An interactive ptx file viewer." OFF)
if(PTXViewer)
  add_subdirectory(GUI)
endif()

option(PTXResectioning "Color the PTX by an external image." OFF)
if(PTXResectioning)
  add_subdirectory(Resectioning)
endif()

# Allow this project to be detected and used as a submodule
set(PTXTools_include_dirs ${PTXTools_include_dirs} ${CMAKE_CURRENT_SOURCE_DIR})
set_property(GLOBAL PROPERTY PTXToolsIncludeDirs ${PTXTools_include_dirs})
set_property(GLOBAL PROPERTY PTXToolsLibs ${PTXTools_libraries})

####### Mass toggle options #######
option(TurnOnAllPTX "Turn on all variables starting with PTX_.")
option(TurnOffAllPTX "Turn off all variables starting with PTX_")

# A function to get all user defined variables with a specified prefix
function (getListOfVarsStartingWith _prefix _varResult)
    get_cmake_property(_vars CACHE_VARIABLES)
    string (REGEX MATCHALL "(^|;)${_prefix}[A-Za-z0-9_]*" _matchedVars "${_vars}")
    set (_resultVars "")
    foreach (_variable ${_matchedVars})
        get_property(_type CACHE "${_variable}" PROPERTY TYPE)
        if (NOT "${_type}" STREQUAL "STATIC")
            list (APPEND _resultVars "${_variable}")
        endif()
    endforeach()

    set (${_varResult} ${_resultVars} PARENT_SCOPE)
endfunction()

if(TurnOnAllPTX AND TurnOffAllPTX)
  message(FATAL_ERROR "You cannot turn targets on and off at the same time!")
endif()

if(TurnOnAllPTX)
    getListOfVarsStartingWith("PTX_" matchedVars)
    foreach (_var IN LISTS matchedVars)
        get_property(currentHelpString CACHE "${_var}" PROPERTY HELPSTRING)
        set(${_var} ON CACHE BOOL "${currentHelpString}" FORCE)
    endforeach()

    get_property(turnOnHelpString CACHE "${TurnOnAllPTX}" PROPERTY HELPSTRING)
    set(TurnOnAllPTX OFF CACHE BOOL "${turnOnHelpString}" FORCE) # Set itself back to off, as this is a one time thing.
endif()

if(TurnOffAllPTX)
    getListOfVarsStartingWith("PTX_" matchedVars)
    foreach (_var IN LISTS matchedVars)
        get_property(currentHelpString CACHE "${_var}" PROPERTY HELPSTRING)
        set(${_var} OFF CACHE BOOL "${currentHelpString}" FORCE)
    endforeach()

    get_property(turnOffHelpString CACHE "${TurnOffAllPTX}" PROPERTY HELPSTRING)
    set(TurnOffAllPTX OFF CACHE BOOL "${turnOffHelpString}" FORCE) # Set itself back to off, as this is a one time thing.
endif()

# Display where this code was used from (if it is used as a submodule, there may be multiple instances of this submodule in the project, only the first of which is used)
option(PTXTools_ShowSubmoduleLocation "Show the path from which PTXTools was used?" OFF)
if(PTXTools_ShowSubmoduleLocation)
  message("PTXTools used from ${CMAKE_CURRENT_SOURCE_DIR}")
endif(PTXTools_ShowSubmoduleLocation)

# Build the tests if requested
option(PTXTools_BuildTests "PTXTools_BuildTests" OFF)
if(PTXTools_BuildTests)
  add_subdirectory(Tests)
endif(PTXTools_BuildTests)
