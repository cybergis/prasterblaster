cmake_minimum_required (VERSION 2.6)
project (pRasterBlaster)

# set release build as default
if (NOT CMAKE_BUILD_TYPE)
		message(STATUS "No build type selected, default to Release")
		set(CMAKE_BUILD_TYPE "Release")
	endif()

# set path to additional CMake modules
set (CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

# set path to bundled GDAL and proj
set (CMAKE_PREFIX_PATH ${CMAKE_SOURCE_DIR}/src/gdal/ ${CMAKE_PREFIX_PATH})

# add src subdirectories to include path
include_directories ("${CMAKE_SOURCE_DIR}" ${CMAKE_SOURCE_DIR}/src/gtest/include/)

# Set CXXFLAGS
set (CMAKE_CXX_FLAGS "-D__PRB_SRC_DIR__=${CMAKE_SOURCE_DIR} ${CMAKE_CXX_FLAGS} -W -Wall -Wextra -Wcast-align -Wpointer-arith -Wsign-compare -Wformat=2 -Wno-format-y2k  -Wmissing-braces -Wparentheses -Wtrigraphs -Wstrict-aliasing=2")

find_package (GDAL)
find_package (Proj)
find_package (TIFF 3.9.4)

if (NOT PROJ_FOUND OR NOT GDAL_FOUND)
  message (SEND_ERROR "Some dependencies were not found!")
endif (NOT PROJ_FOUND OR NOT GDAL_FOUND)

include_directories (${GDAL_INCLUDE_DIR} ${PROJ_INCLUDE_DIR} ${TIFF_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/src/gdal/gdal-1.9.2/frmts/gtiff/libtiff/)

add_library (sptw SHARED src/demos/sptw.cc)
add_library (rasterblaster SHARED src/configuration.cc src/rastercoordtransformer.cc 
  src/reprojection_tools.cc src/rasterchunk.cc)
add_library (prasterblaster SHARED src/demos/prasterblaster-pio.cc)
target_link_libraries (prasterblaster rasterblaster sptw)

add_executable (prasterblasterpio src/demos/prasterblaster-main.cc)
target_link_libraries (prasterblasterpio sptw rasterblaster prasterblaster ${PROJ_LIBRARY} ${GDAL_LIBRARY} ${TIFF_LIBRARY})

add_executable (prasterblaster-simple src/demos/prasterblaster-simple.cc)
target_link_libraries (prasterblaster-simple rasterblaster prasterblaster ${PROJ_LIBRARY} ${GDAL_LIBRARY} ${TIFF_LIBRARY})

add_subdirectory (src/gtest/)

add_executable (tests tests/systemtest.cc tests/check_reprojection_tools.cc tests/rastercompare.cc)

target_link_libraries(tests gtest rasterblaster sptw prasterblaster ${GDAL_LIBRARY} ${PROJ_LIBRARY})

# add a target to generate API documentation with Doxygen
find_package(Doxygen 1.8)
if(DOXYGEN_FOUND)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
  add_custom_target(doc
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen" VERBATIM
    )
endif(DOXYGEN_FOUND)

# add a target to run the 'tests' binary
add_custom_target(check ./tests
                  DEPENDS tests)
