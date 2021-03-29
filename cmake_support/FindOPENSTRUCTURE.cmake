#-------------------------------------------------------------------------------
# Check for OpenStructure Libraries
#
#    OST_ROOT                 Prefix for OpenStructure libraries
#    OST_MIN_VERSION          minimal OpenStructure version required
#  
# When OpenStructure is found, the result is placed in the following variables:
# 
#    OST_LIBRARIES              is set to the library and linker flags used to
#                               link against python
#    OST_VERSION                is set to the version of OpenStructure
#    OST_INCLUDE_DIR            is set to the path that contains base.hh
#    OST_PYMOD_PATH             path to the Pyhton module
#    OST_COMPOUNDS_CHEMLIB_PATH path to the compounds.chemlib
#
# Author: Valerio Mariani, Marco Biasini, Stefan Bienert
#-------------------------------------------------------------------------------

macro(find_OPENSTRUCTURE OST_ROOT HEADER_NAMES PYMOD_NAME)
  if(NOT OPENSTRUCTURE_FIND_COMPONENTS)
    message(FATAL_ERROR "Please specify which modules of OpenStructure you "
            "would like to use after the COMPONENTS keyword.")
  endif()
  list(APPEND OPENSTRUCTURE_FIND_COMPONENTS base geom)
  list(REMOVE_DUPLICATES OPENSTRUCTURE_FIND_COMPONENTS)
  foreach (LIB ${OPENSTRUCTURE_FIND_COMPONENTS})
    set(FOUND_LIB FOUND_LIB-NOTFOUND)
    find_library(FOUND_LIB 
      NAMES ost_${LIB}
      HINTS "${Python_ROOT_DIR}"
      PATH ${OST_ROOT}
      PATH_SUFFIXES lib lib64
      NO_SYSTEM_ENVIRONMENT_PATH NO_DEFAULT_PATH
    )
    if(NOT FOUND_LIB)
      if(OPENSTRUCTURE_FIND_REQUIRED)
        message(FATAL_ERROR "Could not find library ost_${LIB}. Please specify"
                " the location of your OpenStructure installation with"
                " OST_ROOT")
      endif()
    else()
       set(OST_LIBRARIES ${OST_LIBRARIES}  ${FOUND_LIB})
    endif()

  endforeach()

  find_path(OST_INCLUDE_DIR
    NAMES "${HEADER_NAMES}"
    HINTS "${OST_ROOT}/include"
    NO_SYSTEM_ENVIRONMENT_PATH NO_DEFAULT_PATH
    )

  find_path_recursive(OST_PYMOD_PATH
    NAME "${PYMOD_NAME}"
    PATH "${OST_ROOT}"
    PATH_SUFFIXES lib lib64
    )
  if(NOT OST_PYMOD_PATH)
    if(OPENSTRUCTURE_FIND_REQUIRED)
      message(FATAL_ERROR "Could not find Python module of OST. "
                          "Please specify the location of your OST "
                          "installation with OST_ROOT")
    endif(OPENSTRUCTURE_FIND_REQUIRED)
  endif(NOT OST_PYMOD_PATH)

  # compounds.chemlib
  # if no -DCOMPOUND_LIB is given, try to find it in OST_ROOT
  if(COMPOUND_LIB)
    if(NOT EXISTS "${COMPOUND_LIB}")
      message(FATAL_ERROR
        "Could not find compound library at '${COMPOUND_LIB}'")
    endif(NOT EXISTS "${COMPOUND_LIB}")
    if(NOT IS_ABSOLUTE "${COMPOUND_LIB}")
      get_filename_component(COMPOUND_LIB "${COMPOUND_LIB}" ABSOLUTE)
    endif(NOT IS_ABSOLUTE "${COMPOUND_LIB}")
    set(OST_COMPOUNDS_CHEMLIB_PATH "\\\\\\\\\"${COMPOUND_LIB}\\\\\\\\\"")
  else(COMPOUND_LIB)
    set(_clib_name "compounds.chemlib")
    find_path_recursive(OST_COMPOUNDS_CHEMLIB_PATH
      NAME "${_clib_name}"
      PATH "${OST_ROOT}"
      PATH_SUFFIXES share
      )
    if(NOT OST_COMPOUNDS_CHEMLIB_PATH)
      message(WARNING
        "No compound library found (specify with -DCOMPOUND_LIB).")
      set(OST_COMPOUNDS_CHEMLIB_PATH "None")
    else(NOT OST_COMPOUNDS_CHEMLIB_PATH)
      file(TO_NATIVE_PATH "${OST_COMPOUNDS_CHEMLIB_PATH}/${_clib_name}"
        OST_COMPOUNDS_CHEMLIB_PATH)
      set(OST_COMPOUNDS_CHEMLIB_PATH "\\\\\\\\\"${OST_COMPOUNDS_CHEMLIB_PATH}\\\\\\\\\"")
    endif(NOT OST_COMPOUNDS_CHEMLIB_PATH)
  endif(COMPOUND_LIB)
  
  set(OPENSTRUCTURE_FOUND 1)
endmacro(find_OPENSTRUCTURE)

#-------------------------------------------------------------------------------

find_OPENSTRUCTURE("${OST_ROOT}" "ost/config.hh" "ost/__init__.py")

mark_as_advanced(
  OST_LIBRARIES
  OST_INCLUDE_DIR
  OST_VERSION
  OST_PYMOD_PATH
  OST_COMPOUNDS_CHEMLIB_PATH
)

if(OPENSTRUCTURE_FOUND)
   if(NOT OPENSTRUCTURE_FIND_QUIETLY)
   endif()
else(OPENSTRUCTURE_FOUND)
   if(OPENSTRUCTURE_FIND_REQUIRED)
      message(FATAL_ERROR "Could not find OpenStructure")
   endif()
endif()
