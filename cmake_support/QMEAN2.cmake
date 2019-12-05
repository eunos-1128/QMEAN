#-------------------------------------------------------------------------------
# Author: Marco Biasini, Juergen Haas
#
# This file contains a bunch of useful macros to facilitate the build-system
# configuration for the modules.
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# map macro
#
# this function emulates a map/dict data type
#-------------------------------------------------------------------------------

function(map COMMAND MAPNAME)
  set (_KEYS ${MAPNAME}_MAP_KEYS )
  set (_VALUES ${MAPNAME}_MAP_VALUES)
  if(${COMMAND} STREQUAL SET)
    list(REMOVE_AT ARGN 0)
    list(FIND ${_KEYS} ${ARGV2} _MAP_INDEX)
    if(_MAP_INDEX EQUAL -1)
      list(APPEND ${_KEYS} ${ARGV2})
      set(${_KEYS} ${${_KEYS}} PARENT_SCOPE)
      set(${_VALUES}_${ARGV2}  ${ARGN} PARENT_SCOPE)
    else()
      set(${_VALUES}_${ARGV2}  ${ARGN} PARENT_SCOPE)
    endif()
  elseif(${COMMAND} STREQUAL GET)
    list(FIND ${_KEYS} ${ARGV2} _MAP_INDEX)
    if(_MAP_INDEX EQUAL -1)
      MESSAGE(FATAL_ERROR "Unknown key: " ${ARGV2})
    endif()
    set(${ARGV3} ${${_VALUES}_${ARGV2}} PARENT_SCOPE)
  elseif(${COMMAND} STREQUAL KEYS)
    set(${ARGV2} ${${_KEYS}} PARENT_SCOPE)
  elseif(${COMMAND} STREQUAL CREATE)
    set(${_KEYS}  "" PARENT_SCOPE)
  elseif(${COMMAND} STREQUAL LENGTH)
    list(LENGTH ${_KEYS} _L)
    set(${ARGV2} ${_L} PARENT_SCOPE)
  else()
    MESSAGE(FATAL_ERROR "Unknown map command:" ${COMMAND})
  endif()
endfunction()


#-------------------------------------------------------------------------------
# check_architecture
#
# detect architecture based on void pointer size. the output is stored in the
# 3 variables OS_32_BITS, OS_64_BITS and CMAKE_NATIVE_ARCH. The former two are
# set to 0 and 1 accordingly and CMAKE_NATIVE_ARCH is set to the 32 or 64 bit.
#-------------------------------------------------------------------------------
macro(check_architecture)
  include(CheckTypeSize)
  check_type_size(void*  SIZEOF_VOID_PTR)
  if(${SIZEOF_VOID_PTR} MATCHES "^8$")
    set(OS_32_BITS 0)
    set(OS_64_BITS 1)
    set(CMAKE_NATIVE_ARCH 64)
  else()
    set(OS_32_BITS 1)
    set(OS_64_BITS 0)
    set(CMAKE_NATIVE_ARCH 32)
  endif()
endmacro()

#-------------------------------------------------------------------------------
# Synopsis:
#   compile_py_files(module out_dir compiled_files [input_file1 ...])
# Description:
#   Calls pyuic on every input file. The resulting python files are stored in
#   the variable with name compiled_files.
#-------------------------------------------------------------------------------
macro(compile_py_files module out_dir compiled_files_name)
  set(_input_files ${ARGN})
  set(${compiled_files_name})
  foreach(input_file ${_input_files})
    get_filename_component(_out_file ${input_file} NAME_WE)
    get_filename_component(_in_file ${input_file} ABSOLUTE)
    set(_out_file ${out_dir}/${_out_file}.pyc)
    list(APPEND ${compiled_files_name} ${_out_file})
    get_filename_component(_in_name ${input_file} NAME)
    file(MAKE_DIRECTORY  ${out_dir})
    add_custom_command(TARGET ${module}
                       COMMAND ${PYTHON_BINARY} -c "import py_compile;py_compile.compile(\"${_in_file}\",\"${_out_file}\",\"${_in_name}\",doraise=True)"
                       VERBATIM DEPENDS ${input_file}
                       )
  endforeach()
endmacro()

#-------------------------------------------------------------------------------
# this macro has been adapted from
# http://www.cmake.org/Wiki/CMakeMacroParseArguments
#-------------------------------------------------------------------------------
macro(parse_argument_list PREFIX ARG_NAMES OPT_NAMES)
  set(_DEFAULT_ARGS)
  # reset variables
  foreach(${_ARG_NAME} ${ARG_NAMES})
    set(PREFIX_${_ARG_NAME})
  endforeach()
  foreach(_OPT_NAME ${OPT_NAMES})
    set(PREFIX_${_OPT_NAME} FALSE)
  endforeach()
  set(_CURR_ARG_NAME DEF_ARGS)
  set(_CURR_ARG_LIST)
  # loop over parameter list and split by ARG_NAMES
  foreach(_ARG ${ARGN})
    set(_LARG_NAMES ${ARG_NAMES})  
    list(FIND _LARG_NAMES ${_ARG} _IS_ARG_NAME)
    if (_IS_ARG_NAME GREATER -1)
      set(${PREFIX}_${_CURR_ARG_NAME} ${_CURR_ARG_LIST})
      set(_CURR_ARG_NAME "${_ARG}")
      set(_CURR_ARG_LIST)
    else()
    set(_LOPT_NAMES ${OPT_NAMES})  
      list(FIND _LOPT_NAMES ${_ARG} _IS_OPT_NAME)
      if (_IS_OPT_NAME GREATER -1)
        set(${PREFIX}_${_ARG} TRUE)
      else()
        list(APPEND _CURR_ARG_LIST "${_ARG}")
      endif()
    endif()
  endforeach()
  set(${PREFIX}_${_CURR_ARG_NAME} ${_CURR_ARG_LIST})
endmacro()

#-------------------------------------------------------------------------------
# copy_if_different
#
# copies file from source directory to destination directory, but only if its
# content changed.
#-------------------------------------------------------------------------------
macro(copy_if_different FROM_DIR TO_DIR FILES TARGETS TARGET)
  set(ADD_TARGETS "")
  foreach(SRC ${FILES})
      set(SRCFILE ${SRC})
      if("${FROM_DIR}" STREQUAL "" OR "${FROM_DIR}" STREQUAL "./")
          set(FROM ${SRC})
      else()
          set(FROM ${FROM_DIR}/${SRC})
      endif()
      if("${TO_DIR}" STREQUAL "")
          set(TO ${SRCFILE})
      else()
          get_filename_component(TOFILE ${SRC} NAME)      
          set(TO ${TO_DIR}/${TOFILE})
      endif()
      file(MAKE_DIRECTORY  ${TO_DIR})
      #message("copy following stuff: ${FROM} ${TO}")
      add_custom_command(TARGET "${TARGET}" PRE_BUILD
          DEPENDS ${FROM}
          COMMAND ${CMAKE_COMMAND} -E copy_if_different ${FROM} ${TO}
          COMMENT ""
          )
      list(APPEND ADD_TARGETS ${TO})
  endforeach()
  if (TARGETS)
    set(${TARGETS} ${ADD_TARGETS})
  endif()
endmacro()

#-------------------------------------------------------------------------------
# parse_file_list
#
# this macro splits a list of files with IN_DIR statements and fills them into a map
# where the key is the directory name
#-------------------------------------------------------------------------------
macro(parse_file_list FILELIST FILEMAP)
  set(_EXPECT_IN_DIR FALSE)
  map(CREATE ${FILEMAP})
  set(_CURRENT_LIST)
  foreach(_ITEM ${FILELIST})
    if (_ITEM STREQUAL "IN_DIR")
      set(_EXPECT_IN_DIR TRUE)
    else()
      if (_EXPECT_IN_DIR)
        set(_EXPECT_IN_DIR FALSE)
        map(SET ${FILEMAP} ${_ITEM} ${_CURRENT_LIST})
        set(_CURRENT_LIST)
      else()
        list(APPEND _CURRENT_LIST "${_ITEM}")
      endif()
    endif()
  endforeach()
  if(_CURRENT_LIST)
    map(SET ${FILEMAP} "." ${_CURRENT_LIST})
  endif()
endmacro()

#-------------------------------------------------------------------------------
# stage_headers
#-------------------------------------------------------------------------------
macro(stage_headers HEADERS HEADER_INSTALL_DIR TARGET SUB)
  set(FROM_DIR "./")
  # introduce a helper target to make sure the headers are staged before
  # building the library
  string(REPLACE "/" "_" _SUB_NO_SLASH "${SUB}")
  string(REPLACE "${PREFIX}_" "" _TARGET "${TARGET}")
  set(_TARGET_NAME ${_TARGET}_${_SUB_NO_SLASH}_headers)
  set(_SUB ${SUB})
  if (NOT _SUB)
    set(_TARGET_NAME ${_TARGET}_headers)
  endif()
  add_custom_target("${_TARGET_NAME}" COMMENT "")
  set(HEADER_DIR "${HEADER_STAGE_PATH}/${HEADER_INSTALL_DIR}")
  copy_if_different("${FROM_DIR}" "${HEADER_DIR}"
                    "${HEADERS}" ""
                    "${_TARGET_NAME}")
  add_dependencies("${TARGET}" "${_TARGET_NAME}")
  set_target_properties("${_TARGET_NAME}" PROPERTIES
                        EchoString "Staging headers ${TARGET}")
endmacro()


#-------------------------------------------------------------------------------
# Synopsis:
#   module(NAME name SOURCES source1 source2 HEADERS header1 header2 
#          [IN_DIR dir] [header3 header4 [IN_DIR dir]] [DEPENDS_ON dep1 dep2]
#          [HEADER_OUTPUT_DIR dir]
#          [LINK link_cmds])
# Description:
#   Define an OpenStructure module.
#-------------------------------------------------------------------------------
macro(module)
  #-----------------------------------------------------------------------------
  # deal with arguments
  #-----------------------------------------------------------------------------
  set(_ARGS "NAME;SOURCES;HEADERS;DEPENDS_ON;LINK;HEADER_OUTPUT_DIR;PREFIX")
  set(_ARG_PREFIX qmean)  
  parse_argument_list(_ARG "${_ARGS}" "" ${ARGN})  
  if (NOT _ARG_NAME)
    message(FATAL_ERROR 
            "invalid use of module(): a module name must be provided")
  endif()

  if (_ARG_HEADER_OUTPUT_DIR)
    set(_HEADER_OUTPUT_DIR ${_ARG_HEADER_OUTPUT_DIR})
  else()
    if (_ARG_PREFIX)
      set(_HEADER_OUTPUT_DIR "${_ARG_PREFIX}/${_ARG_NAME}")
    else()
      set(_HEADER_OUTPUT_DIR "${_ARG_NAME}")
    endif()
  endif()
  if (_ARG_PREFIX)
    set(_LIB_NAME ${_ARG_PREFIX}_${_ARG_NAME})
  else()
    set(_LIB_NAME ${_ARG_NAME})
  endif()
  string(TOUPPER ${_LIB_NAME} _UPPER_LIB_NAME)  
  #-----------------------------------------------------------------------------
  # create library  
  #-----------------------------------------------------------------------------
  file(MAKE_DIRECTORY ${LIB_STAGE_PATH})
  if (WIN32)
    set(_ABS_FILE_PATTERN "^[A-Z]:/")
  else()
    set(_ABS_FILE_PATTERN "^/")
  endif()
  if (_ARG_SOURCES)
    # when there is at least one source file, we build a library
    set(_ABS_SOURCE_NAMES)
    foreach(_SOURCE ${_ARG_SOURCES})
      if (_SOURCE MATCHES ${_ABS_FILE_PATTERN})
        list(APPEND _ABS_SOURCE_NAMES "${_SOURCE}")
      else()
        list(APPEND _ABS_SOURCE_NAMES "${CMAKE_CURRENT_SOURCE_DIR}/${_SOURCE}")
      endif()
    endforeach()
    add_library(${_LIB_NAME} SHARED ${_ABS_SOURCE_NAMES})
    
    set_target_properties(${_LIB_NAME} 
                          PROPERTIES OUTPUT_NAME ${_LIB_NAME}
                                     PROJECT_LABEL ${_ARG_NAME}
                                     EchoString   ${_ARG_NAME}
                                     MODULE_DEPS "${_ARG_DEPENDS_ON}")
    get_target_property(_DEFS ${_LIB_NAME} COMPILE_DEFINITIONS)
    set_target_properties(${_LIB_NAME} PROPERTIES
                          COMPILE_DEFINITIONS HELMM_MODULE_${_UPPER_LIB_NAME})
    set_target_properties(${_LIB_NAME} PROPERTIES
                          LIBRARY_OUTPUT_DIRECTORY ${LIB_STAGE_PATH}
                          ARCHIVE_OUTPUT_DIRECTORY ${LIB_STAGE_PATH}
                          RUNTIME_OUTPUT_DIRECTORY ${LIB_STAGE_PATH}
                          INSTALL_RPATH "."
                          INSTALL_NAME_DIR "@rpath")
    if (APPLE)
      set_target_properties(${_LIB_NAME} PROPERTIES
                            LINK_FLAGS "-Wl,-rpath,.")
    endif()
    if (WIN32)
      #set_target_properties(${_LIB_NAME} PROPERTIES PREFIX "../")
      install(TARGETS ${_LIB_NAME} ARCHIVE DESTINATION "${LIB_DIR}")
    else()
      install(TARGETS ${_LIB_NAME} LIBRARY DESTINATION "${LIB_DIR}")
    endif()                          
    if (_ARG_LINK)
      target_link_libraries(${_LIB_NAME} ${_ARG_LINK})
    endif()
    foreach(_DEPENDENCY ${_ARG_DEPENDS_ON})
      target_link_libraries(${_LIB_NAME} qmean_${_DEPENDENCY})
    endforeach()
  else()
    add_custom_target("${_LIB_NAME}" ALL)
    set_target_properties("${_LIB_NAME}" PROPERTIES HEADER_ONLY 1 
                          MODULE_DEPS "${_ARG_DEPENDS_ON}")
  endif()
  #-----------------------------------------------------------------------------
  # stage headers  
  #-----------------------------------------------------------------------------
  if (_ARG_HEADERS)
    set(_HEADERS)
    set(_EXPECT_IN_DIR FALSE)
    foreach(_HEADER ${_ARG_HEADERS})
      if (_HEADER STREQUAL "IN_DIR")
        set(_EXPECT_IN_DIR TRUE)
      else()
        if (_EXPECT_IN_DIR)
          set(_EXPECT_IN_DIR FALSE)
          set(_DIR ${_HEADER})
          set(_ABS_HEADER_NAMES)
          set(_HDR_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${_DIR}")
          foreach(_HDR ${_HEADERS})
            list(APPEND _ABS_HEADER_NAMES "${_HDR_SOURCE_DIR}/${_HDR}")
          endforeach()
          install(FILES ${_ABS_HEADER_NAMES} DESTINATION
                  "include/${_HEADER_OUTPUT_DIR}/${_DIR}")
          set(_HDR_STAGE_DIR "${_HEADER_OUTPUT_DIR}/${_DIR}")
          stage_headers("${_ABS_HEADER_NAMES}" "${_HDR_STAGE_DIR}" 
                        "${_LIB_NAME}" "${_DIR}")
          set(_HEADERS)
        else()
          list(APPEND _HEADERS "${_HEADER}")
        endif()
      endif()
    endforeach()
    list(LENGTH _HEADERS _HEADER_LIST_LENGTH)
    if (_HEADER_LIST_LENGTH GREATER 0)    
      set(_ABS_HEADER_NAMES)   
      set(_HDR_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")       
      foreach(_HDR ${_HEADERS})
        list(APPEND _ABS_HEADER_NAMES "${_HDR_SOURCE_DIR}/${_HDR}")
      endforeach()      
      install(FILES ${_ABS_HEADER_NAMES} DESTINATION
              "include/${_HEADER_OUTPUT_DIR}")
      set(_HDR_STAGE_DIR "${_HEADER_OUTPUT_DIR}")
      stage_headers("${_ABS_HEADER_NAMES}" "${_HDR_STAGE_DIR}" 
                    "${_LIB_NAME}" "")
    endif()
  endif()
endmacro()


#-------------------------------------------------------------------------------
# Synopsis
#   executable(NAME exe_name SOURCES source1 source2 LINK link1 link2)
#
# Description:
#  Compile, link and stage a C++ executable
#-------------------------------------------------------------------------------
macro(executable)
  parse_argument_list(_ARG 
                      "NAME;SOURCES;LINK;DEPENDS_ON" "NO_RPATH" ${ARGN})
  if (NOT _ARG_NAME)
    message(FATAL_ERROR "invalid use of executable(): a name must be provided")
  endif()
  add_executable(${_ARG_NAME} ${_ARG_SOURCES})
  if (APPLE AND NOT _ARG_NO_RPATH)
    set_target_properties(${_ARG_NAME} PROPERTIES
                          LINK_FLAGS "-Wl,-rpath,@loader_path/../lib")
  endif()
  if (_ARG_LINK)
    target_link_libraries(${_ARG_NAME} ${_ARG_LINK})
  endif()
  foreach(_DEP ${_ARG_DEPENDS_ON})
    target_link_libraries(${_ARG_NAME} qmean_${_DEP})
  endforeach()
  install(TARGETS ${_ARG_NAME} DESTINATION bin)
endmacro()


#-------------------------------------------------------------------------------
# Synopsis:
#   substitute(IN_FILE in_file OUT_FILE out_file DICT a=b c=d)
#
#-------------------------------------------------------------------------------
macro(substitute)
 parse_argument_list(_ARG 
                     "IN_FILE;OUT_FILE;DICT" "" ${ARGN})
  if (NOT _ARG_IN_FILE)
    message(FATAL_ERROR "invalid use of substitute(): no IN_FILE given")
  endif()
  if (NOT _ARG_OUT_FILE)
    message(FATAL_ERROR "invalid use of substitute(): no OUT_FILE given")
  endif()
  set(_SUBST_DICT -DINPUT_FILE=${_ARG_IN_FILE} -DOUT_FILE=${_ARG_OUT_FILE})
  foreach(_SUBST ${_ARG_DICT})
    list(APPEND _SUBST_DICT -D${_SUBST})
  endforeach()
  add_custom_target(subst_${_ARG_OUT_FILE} ALL COMMAND 
                    ${CMAKE_COMMAND} ${_SUBST_DICT}
                    -P ${CMAKE_SOURCE_DIR}/cmake_support/substitute.cmake)
endmacro()
#-------------------------------------------------------------------------------
# Synopsis:
#   script(NAME script_name INPUT input_name SUBSTITUTE key=val key=val
#         [TARGET target])
#-------------------------------------------------------------------------------
macro(script)
  parse_argument_list(_ARG 
                      "NAME;INPUT;SUBSTITUTE;TARGET" "" ${ARGN})
  if (NOT _ARG_NAME)
    message(FATAL_ERROR "invalid use of executable(): a name must be provided")
  endif()
  set(_INPUT ${_ARG_NAME})
  if (_ARG_INPUT)
    set(_INPUT ${_ARG_INPUT})
  endif()
  if (_ARG_SUBSTITUTE)
    if (NOT _ARG_INPUT)
      message(FATAL_ERROR "script() can only substitute when INPUT is present.")    
    endif()
    substitute(IN_FILE ${_INPUT} OUT_FILE ${_ARG_NAME} 
               DICT ${_ARG_SUBSTITUTE})
  endif()
  install(FILES ${_ARG_NAME} DESTINATION bin 
          PERMISSIONS WORLD_EXECUTE GROUP_EXECUTE OWNER_EXECUTE 
                      WORLD_READ GROUP_READ OWNER_READ)
  copy_if_different("./" "${EXECUTABLE_OUTPUT_PATH}" 
                    "${_ARG_NAME}" "TARGETS" ${_ARG_TARGET})
  add_dependencies(${_ARG_TARGET} subst_${_ARG_NAME})
endmacro()

#-------------------------------------------------------------------------------
# Synopsis:
#   pymod(NAME name CPP source1 source2 PY source source2 [IN_DIR dir] 
#         source3 source4 [IN_DIR dir] [LINK link] [OUTPUT_DIR dir])
#
# Description:
#  Define a python module consisting of C++ type wrappers and/or code written in 
#  Python.
# NAME is the name of
#-------------------------------------------------------------------------------
macro(pymod)
  #-----------------------------------------------------------------------------
  # deal with arguments
  #-----------------------------------------------------------------------------
  set(_ARG_PREFIX qmean)
  parse_argument_list(_ARG 
                      "NAME;CPP;PY;LINK;OUTPUT_DIR;UI;PREFIX" "" ${ARGN})
  if (NOT _ARG_NAME)
    message(FATAL_ERROR "invalid use of pymod(): a name must be provided")
  endif()
  if (ENABLE_STATIC)
    return()
  endif()
  if (_ARG_OUTPUT_DIR)
    set(PYMOD_DIR "${PYTHON_MODULE_PATH}${_ARG_OUTPUT_DIR}")
  else()
    set(PYMOD_DIR "${PYTHON_MODULE_PATH}${_ARG_PREFIX}/${_ARG_NAME}")
  endif()
  if (_ARG_PREFIX)
    set(_LIB_NAME ${_ARG_PREFIX}_${_ARG_NAME})
  else()
    set(_LIB_NAME ${_ARG_NAME})
  endif()
  
  set(PYMOD_STAGE_DIR "${LIB_STAGE_PATH}/${PYMOD_DIR}")
  file(MAKE_DIRECTORY ${PYMOD_STAGE_DIR})
  include_directories(${PYTHON_INCLUDE_PATH})
  #-----------------------------------------------------------------------------
  # compile and link C++ wrappers
  #-----------------------------------------------------------------------------
  if (_ARG_CPP)
    add_library("_${_LIB_NAME}" MODULE ${_ARG_CPP})
    set_target_properties("_${_LIB_NAME}"
                          PROPERTIES ECHO_STRING
                          "Building Python Module ${_ARG_NAME}")
    if (_ARG_PREFIX)
      set(_PARENT_NAME "${_ARG_PREFIX}_${_ARG_NAME}")
    else()
      set(_PARENT_NAME "${_ARG_NAME}")
    endif()
    get_target_property(_CUSTOM_CHECK "${_PARENT_NAME}" HEADER_ONLY)
    if (NOT _CUSTOM_CHECK)
      set(_PARENT_LIB_NAME "${_PARENT_NAME}")
    endif()

    target_link_libraries("_${_LIB_NAME}" ${_PARENT_LIB_NAME} 
                          ${PYTHON_LIBRARIES} ${BOOST_PYTHON_LIBRARIES})

    set_target_properties("_${_LIB_NAME}"
                          PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PYMOD_STAGE_DIR})
    set_target_properties("_${_LIB_NAME}"
                          PROPERTIES LIBRARY_OUTPUT_DIRECTORY_RELEASE ${PYMOD_STAGE_DIR})
    set_target_properties("_${_LIB_NAME}"
                          PROPERTIES LIBRARY_OUTPUT_DIRECTORY_DEBUG ${PYMOD_STAGE_DIR})

    if (NOT ENABLE_STATIC)
      if (_USE_RPATH)
        string(REGEX REPLACE "/[^/]*" "/.." inv_pymod_path "/${PYMOD_DIR}")
        set_target_properties("_${_LIB_NAME}"
                              PROPERTIES INSTALL_RPATH "$ORIGIN${inv_pymod_path}/")
      else()
        set_target_properties("_${_LIB_NAME}"
                              PROPERTIES INSTALL_RPATH "")
      endif()
    endif()
    if (APPLE)
      file(RELATIVE_PATH _REL_PATH "${PYMOD_STAGE_DIR}" "${LIB_STAGE_PATH}")
      set_target_properties("_${_LIB_NAME}" PROPERTIES
                            LINK_FLAGS "-Wl,-rpath,@loader_path/${_REL_PATH}"
                            INSTALL_NAME_DIR "@rpath")
    endif()                          
    if (NOT WIN32)
      set_target_properties("_${_LIB_NAME}"
                          PROPERTIES PREFIX "")
    else ()
      set_target_properties("_${_LIB_NAME}"
                          PROPERTIES SUFFIX ".pyd")
    endif()
    install(TARGETS "_${_LIB_NAME}" LIBRARY DESTINATION
            "${LIB_DIR}/${PYMOD_DIR}")
  else()
    add_custom_target("_${_LIB_NAME}" ALL)
  endif()

  #-----------------------------------------------------------------------------
  # compile python files
  #-----------------------------------------------------------------------------
  if (_ARG_PY)
    parse_file_list("${_ARG_PY}" _PYFILE_MAP)
    map(KEYS _PYFILE_MAP _PYFILE_MAP_KEYS)
    foreach(_DIR ${_PYFILE_MAP_KEYS})
      map(GET _PYFILE_MAP ${_DIR} _PY_FILES)
      set(_ABS_PY_FILES)
      set(_PY_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${_DIR}")
      foreach(_PY ${_PY_FILES})
        list(APPEND _ABS_PY_FILES "${_PY_SOURCE_DIR}/${_PY}")
      endforeach()
      install(FILES ${_ABS_PY_FILES} DESTINATION "${LIB_DIR}/${PYMOD_DIR}/${_DIR}")
      string(REPLACE "/" "_" _DIR_NO_SLASH "${_DIR}")
      set(_PYMOD_TARGET "${_LIB_NAME}_${_DIR_NO_SLASH}_pymod")
      string(REPLACE "_." "" _PYMOD_TARGET "${_PYMOD_TARGET}")
      add_custom_target(${_PYMOD_TARGET} ALL)
      copy_if_different("./" "${PYMOD_STAGE_DIR}/${_DIR}"
                        "${_ABS_PY_FILES}" "TARGETS"
                        "${_PYMOD_TARGET}")
      compile_py_files(_${_LIB_NAME} ${PYMOD_STAGE_DIR}/${_DIR} compiled_files ${_ABS_PY_FILES})
      install(FILES ${compiled_files} DESTINATION "${LIB_DIR}/${PYMOD_DIR}/${_DIR}")
    endforeach()
  endif()  
  get_target_property(_MOD_DEPS "${_PARENT_NAME}" MODULE_DEPS)
  if(_MOD_DEPS)
    foreach(dep ${_MOD_DEPS})
       add_dependencies("_${_LIB_NAME}" "_${dep}")
    endforeach()
  endif()
endmacro()

#-------------------------------------------------------------------------------
# Synopsis:
#   add_unit_test_data_target(DAT target)
#
# Description:
#   Add a data target for a unit test to QM_UNIT_TEST_DATA.
#   DAT - target
#-------------------------------------------------------------------------------
macro(add_unit_test_data_target)
  parse_argument_list(_ARG "DAT" "" ${ARGN})
  if(NOT _ARG_DAT)
    message(FATAL_ERROR
            "invalid use of add_unit_test_data_target(): target missing")
  endif()
  if(DEFINED QM_UNIT_TEST_DATA)
    set(_DATATARGETS "${QM_UNIT_TEST_DATA}")
    set(_NME_FOUND)
    foreach(nme ${QM_UNIT_TEST_DATA})
      if(${nme} STREQUAL ${_ARG_DAT})
        set(_NME_FOUND TRUE)
      endif()
    endforeach()
    if(NOT _NME_FOUND)
      list(APPEND _DATATARGETS "${_ARG_DAT}")
    endif()
  else()
    set(_DATATARGETS "${_ARG_DAT}")
  endif()
  set(QM_UNIT_TEST_DATA "${_DATATARGETS}" CACHE INTERNAL "" FORCE)
endmacro(add_unit_test_data_target)

add_custom_target(check)
if (WIN32)
  set_target_properties(check PROPERTIES EXCLUDE_FROM_ALL "1")
endif()
add_custom_target(check_xml)
add_custom_target(codetest)
add_dependencies(check codetest)


#-------------------------------------------------------------------------------
# get_python_path(VAR):
# 
# Set variable VAR (in parent scope) to an updated PYTHONPATH env. variable 
# which includes OST path. This can then be used to call python scripts
# with commands such as: sh -c "PYTHONPATH=${VAR} ${PYTHON_BINARY} ..."
#-------------------------------------------------------------------------------
macro(get_python_path VAR)
  set(${VAR} $ENV{PYTHONPATH})
  if(${VAR})
    set(${VAR} ":${${VAR}}")
  endif(${VAR})
  set(${VAR} "${LIB_STAGE_PATH}/${PYTHON_MODULE_PATH}${${VAR}}")
  set(${VAR} "${OST_PYMOD_PATH}:${${VAR}}")
endmacro(get_python_path)


#-------------------------------------------------------------------------------
# qmean_unittest
#
# define a unit test
#-------------------------------------------------------------------------------
macro(qmean_unittest)
  set(_ARG_PREFIX qmean)
  parse_argument_list(_ARG "MODULE;SOURCES;LINK;DATA;TARGET;BASE_TARGET" "" ${ARGN})
  set(_SOURCES ${_ARG_SOURCES})
  set(CPP_TESTS)
  set(PY_TESTS)
  if(NOT _ARG_MODULE)
    message(FATAL_ERROR
            "invalid use of qmean_unittest(): module name missing")
  endif()
  set(CMAKE_CURRENT_BINARY_DIR "${CMAKE_BINARY_DIR}/tests/${_ARG_MODULE}")
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  add_custom_target("test_data_${_ARG_MODULE}")
  if(_ARG_TARGET)
    add_dependencies("test_data_${_ARG_MODULE}" "${_ARG_TARGET}")
  endif(_ARG_TARGET)
  add_custom_command(TARGET "test_data_${_ARG_MODULE}" PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR})
  if(_ARG_DATA)
    foreach(dfile ${_ARG_DATA})
      get_filename_component(tdata_path "${dfile}" PATH)
      set(tdata_full_path "${CMAKE_CURRENT_BINARY_DIR}/${tdata_path}")
      if(tdata_path)
        file(MAKE_DIRECTORY "${tdata_full_path}")
        add_custom_command(TARGET "test_data_${_ARG_MODULE}" PRE_BUILD
          COMMAND ${CMAKE_COMMAND} -E make_directory
          "${tdata_full_path}")
      endif()
      get_filename_component(tdata_file "${dfile}" NAME)
      copy_if_different("${CMAKE_CURRENT_SOURCE_DIR}/${tdata_path}"
        "${tdata_full_path}" "${tdata_file}"
        TARGETS "test_data_${_ARG_MODULE}")
      add_unit_test_data_target(DAT "test_data_${_ARG_MODULE}")
    endforeach()
  endif(_ARG_DATA)
  foreach(src ${_SOURCES})
    if(${src} MATCHES "\\.py$")
      list(APPEND PY_TESTS "${src}")
    else()
      list(APPEND CPP_TESTS "${CMAKE_CURRENT_SOURCE_DIR}/${src}")
   endif()
  endforeach()
  set(_SOURCES ${CPP_TESTS})
  set(_test_name "test_suite_${_ARG_MODULE}")
  if(DEFINED CPP_TESTS)
    if(COMPILE_TESTS)
      add_executable(${_test_name} ${_SOURCES})
    else()
      add_executable(${_test_name} EXCLUDE_FROM_ALL ${_SOURCES})
    endif()
    set_target_properties(${_test_name} PROPERTIES
                         RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_target_properties(${_test_name} PROPERTIES
                   RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_CURRENT_BINARY_DIR})
    set_target_properties(${_test_name} PROPERTIES
                 RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_CURRENT_BINARY_DIR})
    add_dependencies(${_test_name} "test_data_${_ARG_MODULE}")
    target_link_libraries(${_test_name} ${BOOST_UNIT_TEST_LIBRARIES}
                        "${_ARG_PREFIX}_${_ARG_MODULE}")
    add_custom_target("${_test_name}_run"
                    COMMAND
           QMEAN_ROOT=${STAGE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/${_test_name}
                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                    COMMENT "running checks for module ${_ARG_MODULE}"
                    DEPENDS ${_test_name})
    if(TARGET "_${_ARG_MODULE}")
      add_dependencies("${_test_name}_run" "_${_ARG_MODULE}")
    endif()
    set(_xml_test_cmd "QMEAN_ROOT=${STAGE_DIR}")
    set(_xml_test_cmd ${_xml_test_cmd} ${CMAKE_CURRENT_BINARY_DIR})
    set(_xml_test_cmd "${_xml_test_cmd}/${_test_name}")
    set(_xml_test_cmd ${_xml_test_cmd} "--log_format=xml" "--log_level=all")
    # XML test outputgets an logical OR to 'echo' so if sth fails, make
    # continues and we get output for all unit tests. Just calling 'echo'
    # giveth $?=0.
    add_custom_target("${_test_name}_run_xml"
                    COMMAND ${_xml_test_cmd} > ${_test_name}_log.xml || echo
                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                    COMMENT "running checks for module ${_ARG_MODULE}"
                    DEPENDS ${_test_name})
    if(TARGET "_${_ARG_MODULE}")
      add_dependencies("${_test_name}_run_xml" "_${_ARG_MODULE}")
    endif()
    add_test("${_test_name}" ${CMAKE_CURRENT_BINARY_DIR}/${_test_name} )
    add_dependencies(check_xml "${_test_name}_run_xml")
    if(_ARG_BASE_TARGET)
      add_dependencies("${_ARG_BASE_TARGET}" "${_test_name}_run")
    else()
      add_dependencies(codetest "${_test_name}_run")
    endif()

    if(_ARG_LINK)
      target_link_libraries("${_test_name}" ${_ARG_LINK})
    endif()

    set_target_properties(${_test_name}
                          PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                          "${CMAKE_CURRENT_BINARY_DIR}")
  endif()

  foreach(py_test ${PY_TESTS})
    set(py_twp "${CMAKE_CURRENT_SOURCE_DIR}/${py_test}")
    if(NOT EXISTS "${py_twp}")
      message(FATAL_ERROR "Python test script does not exist: ${py_twp}")
    endif()

    get_python_path(python_path)

    set (PY_TESTS_CMD "PYTHONPATH=${python_path} QMEAN_PYTHON_BINARY=${PYTHON_BINARY}  ${PYTHON_BINARY}")

    add_custom_target("${py_test}_run"
                sh -c "${PY_TESTS_CMD} ${py_twp}"
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
              COMMENT "running checks ${py_test}" VERBATIM)
    add_dependencies("${py_test}_run" "test_data_${_ARG_MODULE}")
    if(TARGET "_${_ARG_MODULE}")
      add_dependencies("${py_test}_run" "_${_ARG_MODULE}")
    endif()
    # XML test output gets an logical OR to 'echo' so if sth fails, make
    # continues and we get output for all unit tests. Just calling 'echo'
    # gives $?=0.
    add_custom_target("${py_test}_run_xml"
                sh -c "${PY_TESTS_CMD} ${py_twp} xml || echo"
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
              COMMENT "running checks ${py_test}" VERBATIM)
    add_dependencies("${py_test}_run_xml" "test_data_${_ARG_MODULE}")
    if(TARGET "_${_ARG_MODULE}")
      add_dependencies("${py_test}_run_xml" "_${_ARG_MODULE}")
    endif()
    add_dependencies(check_xml "${py_test}_run_xml")
    if(_ARG_BASE_TARGET)
      add_dependencies("${_ARG_BASE_TARGET}" "${py_test}_run")
    else()
      add_dependencies(codetest "${py_test}_run")
    endif()
  endforeach()
endmacro(qmean_unittest)

#-------------------------------------------------------------------------------
# make sure the previously detected Python interpreter has the given module
#-------------------------------------------------------------------------------
macro(qmean_find_python_module MODULE)
  if (NOT PYTHON_MODULE_${MODULE})
    execute_process(COMMAND ${PYTHON_BINARY} -c "import ${MODULE}"
                    OUTPUT_QUIET ERROR_QUIET
                    RESULT_VARIABLE _IMPORT_ERROR)
    if (_IMPORT_ERROR)
      message(FATAL_ERROR "Could not find python module ${MODULE}. Please install it")
    else()
      message(STATUS "Found python module ${MODULE}")
      set("PYTHON_MODULE_${MODULE}" FOUND CACHE STRING "" FORCE)
    endif()
  endif()
endmacro()


#-------------------------------------------------------------------------------
# this macro tries to detect a very common problem during configuration stage:
# it makes sure that boost python is linked against the same version of python 
# that we are linking against. 
#-------------------------------------------------------------------------------
macro(qmean_match_boost_python_version)
  include(CheckCXXSourceRuns)
  # this variable may either be a simple library path or list that contains
  # different libraries for different build-options. For example:
  # optimized;<lib1>;debug;<lib2>
  set(_BOOST_PYTHON_LIBRARY ${Boost_PYTHON_LIBRARY})
  list(LENGTH _BOOST_PYTHON_LIBRARY _BP_LENGTH)
  if (_BP_LENGTH GREATER 1)
    list(FIND _BOOST_PYTHON_LIBRARY optimized _OPTIMIZED_INDEX)
    if (_OPTIMIZED_INDEX EQUAL -1)
      message(FATAL_ERROR 
              "Error while trying to get path of boost python library")
    endif()
    math(EXPR _LIB_INDEX 1+${_OPTIMIZED_INDEX})
    list(GET _BOOST_PYTHON_LIBRARY ${_LIB_INDEX} _BP_LIB_PATH)
    set(_BOOST_PYTHON_LIBRARY ${_BP_LIB_PATH})
  endif()
  set(CMAKE_REQUIRED_FLAGS "-I${PYTHON_INCLUDE_PATH} -I${Boost_INCLUDE_DIR} ${PYTHON_LIBRARIES} ${_BOOST_PYTHON_LIBRARY}")
  check_cxx_source_runs(
"#include <boost/python.hpp>

using namespace boost::python;
                           
int main( int argc, char ** argv ) {
  try {
    Py_Initialize();

    object main_module((
      handle<>(borrowed(PyImport_AddModule(\"__main__\")))));

    object main_namespace = main_module.attr(\"__dict__\");

    handle<> ignored(( PyRun_String( \"print 'Hello, World'\",
                                     Py_file_input,
                                     main_namespace.ptr(),
                                     main_namespace.ptr() ) ));
  } catch( error_already_set ) {
    PyErr_Print();
    return 1;
  }   
  return 0;
}" PYTHON_VERSION_FOUND_MATCHES_PYTHON_VERSION_BOOST_WAS_COMPILED_WITH)

  if(NOT PYTHON_VERSION_FOUND_MATCHES_PYTHON_VERSION_BOOST_WAS_COMPILED_WITH)
     message(FATAL_ERROR "Python versions mismatch!\n"
             "Make sure the python version for boost match the one you "
             "specified during configuration\n")
  endif()
endmacro()

macro(qmean_action ACTION TARGET)
  copy_if_different("./" "${STAGE_DIR}/libexec" 
                    "${ACTION}" "TARGETS" ${TARGET})
  install(FILES "${ACTION}" DESTINATION "libexec")
endmacro()

#-------------------------------------------------------------------------------
# Synopsis:
#   find_path_recursive(VARIABLE
#                       NAME file
#                       PATH path
#                       [PATH_SUFFIXES dir1 dir2 ...])
# Description:
#   Find a path to a file in a directory tree. The result is stored in the
#   given variable. Its 'VARIABLE-NOTFOUND' if the file could not be found.
#   PATH defines were to search, PATH_SUFFIXES sub-directories in PATH to limit
#   the search.
#-------------------------------------------------------------------------------
macro(find_path_recursive VARIABLE)
  set(_ARGS "NAME;PATH;PATH_SUFFIXES")
  parse_argument_list(_ARG "${_ARGS}" "" ${ARGN})
  if(NOT _ARG_NAME)
    message(FATAL_ERROR "Wrong use of find_path_recursive, we need a file NAME "
                        "to look for.")
  endif(NOT _ARG_NAME)
  if(NOT _ARG_PATH)
    message(FATAL_ERROR "Wrong use of find_path_recursive, we need a PATH to "
                        "look in.")
  endif(NOT _ARG_PATH)
  # get first level of dirs (needed to deal with PATH_SUFFIXES if present)
  if(_ARG_PATH_SUFFIXES)
    set(_fst_subs)
    foreach(_psuf ${_ARG_PATH_SUFFIXES})
      file(TO_NATIVE_PATH "${_ARG_PATH}/${_psuf}" _psuf_path)
      if(EXISTS ${_psuf_path})
        list(APPEND _fst_subs ${_psuf_path})
      endif(EXISTS ${_psuf_path})
    endforeach()
  else(_ARG_PATH_SUFFIXES)
    file(GLOB _fst_subs ${_ARG_PATH}/*)
  endif(_ARG_PATH_SUFFIXES)
  set(${VARIABLE} "${VARIABLE}-NOTFOUND")
  # iterate, top-down, until first hit
  while(_fst_subs)
    foreach(_fst ${_fst_subs})
      # check if exists and exit
      file(TO_NATIVE_PATH "${_fst}/${_ARG_NAME}" _modpath)
      if(EXISTS ${_modpath})
        # set var and exit
        set(${VARIABLE} ${_fst})
        set(_fst_subs)
        break()
      endif(EXISTS ${_modpath})
    endforeach(_fst ${_fst_subs})
    set(_tmp_dlist)
    foreach(_fst ${_fst_subs})
      # list subdirs
      file(GLOB _snd_subs ${_fst}/*)
      foreach(_snd ${_snd_subs})
        if(IS_DIRECTORY ${_snd})
          list(APPEND _tmp_dlist ${_snd})
        endif()
      endforeach()
    endforeach(_fst ${_fst_subs})
    # swap lists
    set(_fst_subs ${_tmp_dlist})
  endwhile(_fst_subs)
endmacro(find_path_recursive)
