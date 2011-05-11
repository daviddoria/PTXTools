cmake_minimum_required(VERSION 2.6)

PROJECT(PTXReader)

FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

SET( INSTALL_DIR ${CMAKE_INSTALL_PREFIX} )

if(CMAKE_BUILD_TYPE STREQUAL "Release")
  message("Compiled in Release mode")
else(CMAKE_BUILD_TYPE STREQUAL "Release")
  list(APPEND DEFINITIONS "DEBUGMODE")
  message("Compiled in Debug mode")
endif(CMAKE_BUILD_TYPE STREQUAL "Release")

#set(SourceFiles PTXPixel.cxx PTXImage.cxx)

add_library(PTX PTXPixel.cxx PTXImage.cxx Helpers.cpp)
set(PTXLibraries PTX ${ITK_LIBRARIES} vtkHybrid vtkIO)

#############
ADD_EXECUTABLE(CoordinateImages CoordinateImages.cxx)
TARGET_LINK_LIBRARIES(CoordinateImages ${PTXLibraries})
INSTALL( TARGETS CoordinateImages
  RUNTIME DESTINATION ${INSTALL_DIR} )

ADD_EXECUTABLE(AppendRight AppendRight.cxx)
TARGET_LINK_LIBRARIES(AppendRight ${PTXLibraries})
INSTALL( TARGETS AppendRight
  RUNTIME DESTINATION ${INSTALL_DIR} )

ADD_EXECUTABLE(CropPTX CropPTX.cxx)
TARGET_LINK_LIBRARIES(CropPTX ${PTXLibraries})
INSTALL( TARGETS CropPTX
  RUNTIME DESTINATION ${INSTALL_DIR} )

# ADD_EXECUTABLE(DerivativeImages DerivativeImages.cxx)
# TARGET_LINK_LIBRARIES(DerivativeImages ${PTXLibraries})
# INSTALL( TARGETS DerivativeImages
#   RUNTIME DESTINATION ${INSTALL_DIR} )

ADD_EXECUTABLE(LaplacianImages LaplacianImages.cxx)
TARGET_LINK_LIBRARIES(LaplacianImages ${PTXLibraries})
INSTALL( TARGETS LaplacianImages
  RUNTIME DESTINATION ${INSTALL_DIR} )

# ADD_EXECUTABLE(ExtractHoleBoundaryPoints ExtractHoleBoundaryPoints.cxx)
# TARGET_LINK_LIBRARIES(ExtractHoleBoundaryPoints ${PTXLibraries})
# INSTALL( TARGETS ExtractHoleBoundaryPoints
#   RUNTIME DESTINATION ${INSTALL_DIR} )


ADD_EXECUTABLE(PTXtoRGB PTXtoRGB.cxx)
TARGET_LINK_LIBRARIES(PTXtoRGB ${PTXLibraries})
INSTALL( TARGETS PTXtoRGB
  RUNTIME DESTINATION ${INSTALL_DIR} )


ADD_EXECUTABLE(PTXtoPointCloud PTXtoPointCloud.cxx)
TARGET_LINK_LIBRARIES(PTXtoPointCloud ${PTXLibraries})
INSTALL( TARGETS PTXtoPointCloud
  RUNTIME DESTINATION ${INSTALL_DIR} )


# ADD_EXECUTABLE(PTXtoDepthImage PTXtoDepthImage.cxx)
# TARGET_LINK_LIBRARIES(PTXtoDepthImage ${PTXLibraries})
# INSTALL( TARGETS PTXtoDepthImage
#   RUNTIME DESTINATION ${INSTALL_DIR} )


# ADD_EXECUTABLE(PTXtoDepthThresholdMask PTXtoDepthThresholdMask.cxx)
# TARGET_LINK_LIBRARIES(PTXtoDepthThresholdMask ${PTXLibraries})
# INSTALL( TARGETS PTXtoDepthThresholdMask
#   RUNTIME DESTINATION ${INSTALL_DIR} )


# ADD_EXECUTABLE(ReplaceDepthImage ReplaceDepthImage.cxx)
# TARGET_LINK_LIBRARIES(ReplaceDepthImage ${PTXLibraries})
# INSTALL( TARGETS ReplaceDepthImage
#   RUNTIME DESTINATION ${INSTALL_DIR} )


# ADD_EXECUTABLE(ReplaceRGBD ReplaceRGBD.cxx)
# TARGET_LINK_LIBRARIES(ReplaceRGBD ${PTXLibraries})
# INSTALL( TARGETS ReplaceRGBD
#   RUNTIME DESTINATION ${INSTALL_DIR} )

ADD_EXECUTABLE(ReplaceRGB ReplaceRGB.cxx)
TARGET_LINK_LIBRARIES(ReplaceRGB ${PTXLibraries})
INSTALL( TARGETS ReplaceRGB
  RUNTIME DESTINATION ${INSTALL_DIR} )

ADD_EXECUTABLE(ReplaceXYZ ReplaceXYZ.cxx)
TARGET_LINK_LIBRARIES(ReplaceXYZ ${PTXLibraries})
INSTALL( TARGETS ReplaceXYZ
  RUNTIME DESTINATION ${INSTALL_DIR} )

ADD_EXECUTABLE(PTXtoEverything PTXtoEverything.cxx)
TARGET_LINK_LIBRARIES(PTXtoEverything ${PTXLibraries})
INSTALL( TARGETS PTXtoEverything
  RUNTIME DESTINATION ${INSTALL_DIR} )


ADD_EXECUTABLE(PTXtoValidityMask PTXtoValidityMask.cxx)
TARGET_LINK_LIBRARIES(PTXtoValidityMask ${PTXLibraries})
INSTALL( TARGETS PTXtoValidityMask
  RUNTIME DESTINATION ${INSTALL_DIR} )

# 
# ADD_EXECUTABLE(WeightedLaplacian WeightedLaplacian.cxx)
# TARGET_LINK_LIBRARIES(WeightedLaplacian ${PTXLibraries})
# INSTALL( TARGETS WeightedLaplacian
#   RUNTIME DESTINATION ${INSTALL_DIR} )
# 
# 
# ADD_EXECUTABLE(Downsample Downsample.cxx)
# TARGET_LINK_LIBRARIES(Downsample ${PTXLibraries})
# INSTALL( TARGETS Downsample
#   RUNTIME DESTINATION ${INSTALL_DIR} )

# 
# ADD_EXECUTABLE(TestApplyMask TestApplyMask.cxx)
# TARGET_LINK_LIBRARIES(TestApplyMask ${PTXLibraries})
# INSTALL( TARGETS TestApplyMask
#   RUNTIME DESTINATION ${INSTALL_DIR} )
# 
# 
# ADD_EXECUTABLE(OrthogonalProjection OrthogonalProjection.cxx)
# TARGET_LINK_LIBRARIES(OrthogonalProjection ${PTXLibraries})
# set_target_properties(OrthogonalProjection PROPERTIES COMPILE_DEFINITIONS "${DEFINITIONS}")
# INSTALL( TARGETS OrthogonalProjection
#   RUNTIME DESTINATION ${INSTALL_DIR} )

