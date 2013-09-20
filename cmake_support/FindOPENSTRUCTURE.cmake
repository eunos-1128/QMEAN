#-------------------------------------------------------------------------------
# Check for OpenStructure Libraries
#
#    OST_ROOT                 Prefix for OpenStructure libraries
#    OST_MIN_VERSION          minimal OPenStructure version required
#  
# When OpenStructure is found, the result is placed in the following variables:
# 
#    OST_LIBRARIES           is set to the library and linker flags used to
#                            link against python
#    OST_VERSION             is set to the version of OpenStructure
#    OST_INCLUDE_DIR         is set to the path that contains base.hh
#
# Author: Valerio Mariani, Marco Biasini
#-------------------------------------------------------------------------------

macro(find_OPENSTRUCTURE OST_ROOT NAMES HEADER_NAMES)
  if (NOT OPENSTRUCTURE_FIND_COMPONENTS)
    message(FATAL_ERROR "Please specify which modules of OpenStructure you "
            "would like to use after the COMPONENTS keyword.")
  endif()
  list(APPEND OPENSTRUCTURE_FIND_COMPONENTS base geom)
  list(REMOVE_DUPLICATES OPENSTRUCTURE_FIND_COMPONENTS)
  foreach (LIB ${OPENSTRUCTURE_FIND_COMPONENTS})
    set(FOUND_LIB FOUND_LIB-NOTFOUND)
    find_library(FOUND_LIB 
      NAMES ost_${LIB}
      HINTS "${PYTHON_ROOT}"
      PATH ${OST_ROOT}
      PATH_SUFFIXES lib lib64
      NO_SYSTEM_ENVIRONMENT_PATH NO_DEFAULT_PATH
    )
    if (NOT FOUND_LIB)
      if (OPENSTRUCTURE_FIND_REQUIRED)
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
  set(OPENSTRUCTURE_FOUND 1)
endmacro()

#-------------------------------------------------------------------------------
SET(LIBNAMES
ost_base
ost_gfx
ost_info
ost_qa
ost_conop
ost_gui
ost_io
ost_seq
ost_db
ost_img
ost_mol
ost_seq
ost_seq_alg
ost_geom
ost_img_alg	
ost_mol_alg
)

find_OPENSTRUCTURE("${OST_ROOT}" "${LIBNAMES}" "ost/config.hh")

mark_as_advanced(
  OST_LIBRARIES
  OST_INCLUDE_DIR
  OST_VERSION
)

if (OPENSTRUCTURE_FOUND)
   if (NOT OPENSTRUCTURE_FIND_QUIETLY)
   endif ()
else (OPENSTRUCTURE_FOUND)
   if (OPENSTRUCTURE_FIND_REQUIRED)
      message(FATAL_ERROR "Could not find OpenStructure")
   endif () 
endif ()