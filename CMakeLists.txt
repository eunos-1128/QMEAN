#-------------------------------------------------------------------------------
# Author: Marco Biasini
#-------------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.12.1 FATAL_ERROR)

# Set CMake policies
# Behaviour of target_link_libraries, always link by full path. CMP0060 can
# disable this to stay compatible with old projects. We want the new behaviour
# since before we needed to provide the full path ourselves, that is: we
# manually forced CMake to the new scheme.
cmake_policy(SET CMP0060 NEW)

set (QMEAN_VERSION_MAJOR 4)
set (QMEAN_VERSION_MINOR 1)
set (QMEAN_VERSION_PATCH 0)
set (QMEAN_VERSION_STRING ${QMEAN_VERSION_MAJOR}.${QMEAN_VERSION_MINOR}.${QMEAN_VERSION_PATCH} )

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake_support)

include(QMEAN2)

project(QMEAN2 CXX C)

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

setup_compiler_flags()

# Python needed before Boost
find_package(Python 3.6 REQUIRED COMPONENTS Interpreter Development)
# Split version string
string(REPLACE "." ";" _python_version_list ${Python_VERSION})
list(GET _python_version_list 0 Python_VERSION_MAJOR)
list(GET _python_version_list 1 Python_VERSION_MINOR)
# where Python modules live
set(PYTHON_MODULE_PATH "python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages/")

setup_boost()

find_package(OPENSTRUCTURE 2.1.0 REQUIRED
             COMPONENTS mol seq seq_alg mol_alg conop db)

include_directories(${Boost_INCLUDE_DIRS}
                    ${OST_INCLUDE_DIR})
                    
add_subdirectory(data)
add_subdirectory(src)
add_subdirectory(pymod)
add_subdirectory(doc)

set (FILES_TO_BE_REMOVED stage CMakeFiles)
set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
                         "${FILES_TO_BE_REMOVED}")

message(STATUS "QMEAN will be built with the following options:\n"
        "   OpenStructure (-DOST_ROOT) : ${OST_ROOT}\n"
        "   Optimized     (-DOPTIMIZE) : ${_OPT}\n")
