# Allow headers in tests to be included like
# #include "PTXPixel.h" rather than needing
# #include "PTXTools/PTXPixel.h"
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/..)

add_executable(TestPTXTools TestPTXTools.cpp)
target_link_libraries(TestPTXTools PTX)

add_executable(TestRealPTX TestRealPTX.cpp)
target_link_libraries(TestRealPTX PTX)
