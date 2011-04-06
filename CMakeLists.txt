cmake_minimum_required(VERSION 2.6)

PROJECT(PTXReader)

FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

#set(SourceFiles PTXPixel.cxx PTXImage.cxx)

add_library(PTX PTXPixel.cxx PTXImage.cxx)
set(PTXLibraries PTX ${ITK_LIBRARIES} vtkHybrid vtkIO)

ADD_EXECUTABLE(PTXtoRGB PTXtoRGB.cxx)
TARGET_LINK_LIBRARIES(PTXtoRGB ${PTXLibraries})

ADD_EXECUTABLE(PTXtoPointCloud PTXtoPointCloud.cxx)
TARGET_LINK_LIBRARIES(PTXtoPointCloud ${PTXLibraries})

ADD_EXECUTABLE(PTXtoDepthImage PTXtoDepthImage.cxx)
TARGET_LINK_LIBRARIES(PTXtoDepthImage ${PTXLibraries})

ADD_EXECUTABLE(PTXtoDepthThresholdMask PTXtoDepthThresholdMask.cxx)
TARGET_LINK_LIBRARIES(PTXtoDepthThresholdMask ${PTXLibraries})

ADD_EXECUTABLE(ReplaceDepthImage ReplaceDepthImage.cxx)
TARGET_LINK_LIBRARIES(ReplaceDepthImage ${PTXLibraries})

ADD_EXECUTABLE(ReplaceRGBD ReplaceRGBD.cxx)
TARGET_LINK_LIBRARIES(ReplaceRGBD ${PTXLibraries})

ADD_EXECUTABLE(PTXtoEverything PTXtoEverything.cxx)
TARGET_LINK_LIBRARIES(PTXtoEverything ${PTXLibraries})

ADD_EXECUTABLE(PTXtoValidityMask PTXtoValidityMask.cxx)
TARGET_LINK_LIBRARIES(PTXtoValidityMask ${PTXLibraries})

ADD_EXECUTABLE(WeightedLaplacian WeightedLaplacian.cxx)
TARGET_LINK_LIBRARIES(WeightedLaplacian ${PTXLibraries})

ADD_EXECUTABLE(Downsample Downsample.cxx)
TARGET_LINK_LIBRARIES(Downsample ${PTXLibraries})

ADD_EXECUTABLE(TestApplyMask TestApplyMask.cxx)
TARGET_LINK_LIBRARIES(TestApplyMask ${PTXLibraries})

ADD_EXECUTABLE(OrthogonalProjection OrthogonalProjection.cxx)
TARGET_LINK_LIBRARIES(OrthogonalProjection ${PTXLibraries})