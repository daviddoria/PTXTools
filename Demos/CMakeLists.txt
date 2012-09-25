# Allow headers in tests to be included like
# #include "PTXPixel.h" rather than needing
# #include "PTXTools/PTXPixel.h"
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/..)

add_executable(DemoPTXCoordinateSystem DemoPTXCoordinateSystem.cpp)
target_link_libraries(DemoPTXCoordinateSystem PTX)
