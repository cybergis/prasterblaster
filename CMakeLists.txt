cmake_minimum_required (VERSION 2.8)
project (pRasterBlaster)

# set path to additional CMake modules
set (CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

# set path to bundled GDAL and proj
set (CMAKE_PREFIX_PATH ${CMAKE_SOURCE_DIR}/src/gdal/ ${CMAKE_PREFIX_PATH})

# add src subdirectories to include path
include_directories ("${CMAKE_SOURCE_DIR}")

# Set CXXFLAGS
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -Wextra -Wcast-align -Wpointer-arith -Wsign-compare -Wformat=2 -Wno-format-y2k  -Wmissing-braces -Wparentheses -Wtrigraphs -Wstrict-aliasing=2")

find_package (GDAL)
find_package (Proj)

if (NOT PROJ_FOUND OR NOT GDAL_FOUND)
  message (SEND_ERROR "Some dependencies were not found!")
endif (NOT PROJ_FOUND OR NOT GDAL_FOUND)

include_directories (${GDAL_INCLUDE_DIR} ${PROJ_INCLUDE_DIR})

add_library (sptw src/demos/sptw.cc)
add_library (rasterblaster src/configuration.cc src/quadtree.cc src/rastercoordtransformer.cc 
  src/reprojection_tools.cc src/rasterchunk.cc)

add_executable (prasterblaster src/demos/prasterblaster-pio.cc src/demos/prasterblaster-main.cc)
target_link_libraries (prasterblaster sptw rasterblaster ${PROJ_LIBRARY} ${GDAL_LIBRARY})

add_subdirectory (src/gtest/)

add_executable (tests tests/check_gdal.cc tests/check-rastercoordtransformer.cc tests/check-reprojection_tools.cc)
add_executable (systemtest tests/systemtest.cc src/demos/prasterblaster-pio.cc)

target_link_libraries(tests gtest gtest_main ${GDAL_LIBRARY} ${PROJ_LIBRARY})
target_link_libraries(systemtest rasterblaster sptw gtest ${GDAL_LIBRARY} ${PROJ_LIBRARY})

