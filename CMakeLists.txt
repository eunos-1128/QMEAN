#-------------------------------------------------------------------------------
# Author: Marco Biasini
#-------------------------------------------------------------------------------
cmake_minimum_required(VERSION 2.6.4 FATAL_ERROR)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake_support)

include(QMEAN2)

project(QMEAN2)

if (PREFIX)
  set(CMAKE_INSTALL_PREFIX ${PREFIX})
endif()

if (SYS_ROOT)
  set(CMAKE_PREFIX_PATH ${SYS_ROOT})
endif()

if (OPTIMIZE)
  set(CMAKE_BUILD_TYPE Release)
  set(_OPT ON)
else()
   set(CMAKE_BUILD_TYPE Debug)
    set(_OPT OFF)
endif()

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif()

set(STAGE_DIR "${CMAKE_BINARY_DIR}/stage")
set(HEADER_STAGE_PATH ${STAGE_DIR}/include)
set(SHARED_DATA_PATH ${STAGE_DIR}/share/qmean/)
include_directories("${HEADER_STAGE_PATH}")
link_directories(${LIB_STAGE_PATH})

if (UNIX AND NOT APPLE)
  check_architecture()
endif()
set (ARCH ${CMAKE_NATIVE_ARCH})
if ("${ARCH}" MATCHES "64")
  set(LIB_DIR lib64)
  set(LIB_STAGE_PATH "${STAGE_DIR}/lib64")
else()
  set(LIB_DIR lib)
  set(LIB_STAGE_PATH "${STAGE_DIR}/lib")
endif()

file(MAKE_DIRECTORY ${STAGE_DIR} ${HEADER_STAGE_PATH} ${LIB_STAGE_PATH})

find_package(Python 2.7.0 REQUIRED)

find_package(OPENSTRUCTURE 1.7 REQUIRED 
             COMPONENTS mol seq seq_alg mol_alg conop)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY)
set(_BOOST_MIN_VERSION 1.53)
find_package(Boost ${_BOOST_MIN_VERSION} 
             COMPONENTS filesystem iostreams system REQUIRED)
set(BOOST_LIBRARIES ${Boost_LIBRARIES})
set(Boost_LIBRARIES)
find_package(Boost ${_BOOST_MIN_VERSION} 
             COMPONENTS unit_test_framework REQUIRED)
set(BOOST_UNIT_TEST_LIBRARIES ${Boost_LIBRARIES})
set(Boost_LIBRARIES)
find_package(Boost ${_BOOST_MIN_VERSION} COMPONENTS python REQUIRED)
set(BOOST_PYTHON_LIBRARIES ${Boost_LIBRARIES})
set(Boost_LIBRARIES)

find_package(Python 2.7.0 REQUIRED)
include_directories(${Boost_INCLUDE_DIRS} 
                    ${OST_INCLUDE_DIR})
                    
add_subdirectory(data)
add_subdirectory(src)
add_subdirectory(pymod)

set (FILES_TO_BE_REMOVED stage CMakeFiles)
set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
                         "${FILES_TO_BE_REMOVED}")

message(STATUS "QMEAN will be built with the following options:\n"
        "   OpenStructure (-DOST_ROOT) : ${OST_ROOT}\n"
        "   Optimized     (-DOPTIMIZE) : ${_OPT}\n")
