set(SUBST_DICT BUILD_TYPE=${CMAKE_BUILD_TYPE} 
               PYTHON_BINARY=${PYTHON_BINARY} LIBDIR=${LIB_DIR}
               OST_ROOT=${OST_ROOT})

add_custom_target(qmean_scripts ALL)
script(NAME qmean INPUT qmean.in SUBSTITUTE ${SUBST_DICT}
       TARGET qmean_scripts)

set(INIT_SCRIPTS 
  qmean_script.py
)

copy_if_different("./" "${LIB_EXEC}" 
                  "${INIT_SCRIPTS}" "qmean init script" qmean_scripts)

install(FILES ${INIT_SCRIPTS} DESTINATION "libexec/qmean/")