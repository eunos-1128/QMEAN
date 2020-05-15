set(QMEAN_DATA_DIR "${CMAKE_SOURCE_DIR}/data/qmean")

set(POTENTIAL_DIR "${QMEAN_DATA_DIR}/potentials")
set(POTENTIAL_FILES
${POTENTIAL_DIR}/soluble_global_potentials.dat 
${POTENTIAL_DIR}/soluble_local_potentials.dat 
${POTENTIAL_DIR}/membrane_local_potentials.dat
${POTENTIAL_DIR}/membrane_global_potentials.dat
)

set(REFERENCE_TABLE_DIR "${QMEAN_DATA_DIR}/reference_values")
set(REFERENCE_TABLE_FILES
${REFERENCE_TABLE_DIR}/reference_tab_promod_scorer.txt
)

set(LINEAR_SCORER_DIR "${QMEAN_DATA_DIR}/scorer")
set(LINEAR_SCORER_FILES
${LINEAR_SCORER_DIR}/promod_local_scorer.dat
${LINEAR_SCORER_DIR}/promod_global_scorer.dat
${LINEAR_SCORER_DIR}/membrane_local_scorer.dat
)

set(LOCAL_NN_SCORER_DIR "${QMEAN_DATA_DIR}/scorer/local_nn_scorer")
set(LOCAL_NN_SCORER_FILES
${LOCAL_NN_SCORER_DIR}/feature_groups.json
${LOCAL_NN_SCORER_DIR}/nn_0.dat
${LOCAL_NN_SCORER_DIR}/nn_1.dat
${LOCAL_NN_SCORER_DIR}/nn_2.dat
${LOCAL_NN_SCORER_DIR}/nn_3.dat
${LOCAL_NN_SCORER_DIR}/nn_4.dat
${LOCAL_NN_SCORER_DIR}/nn_5.dat
${LOCAL_NN_SCORER_DIR}/nn_6.dat
${LOCAL_NN_SCORER_DIR}/nn_7.dat
${LOCAL_NN_SCORER_DIR}/nn_8.dat
${LOCAL_NN_SCORER_DIR}/nn_9.dat
${LOCAL_NN_SCORER_DIR}/nn_10.dat
${LOCAL_NN_SCORER_DIR}/nn_11.dat
${LOCAL_NN_SCORER_DIR}/nn_12.dat
${LOCAL_NN_SCORER_DIR}/nn_13.dat
${LOCAL_NN_SCORER_DIR}/nn_14.dat
${LOCAL_NN_SCORER_DIR}/nn_15.dat
)

set(GMQE_SCORER_DIR "${QMEAN_DATA_DIR}/scorer/gmqe_scorer")
set(GMQE_SCORER_FILES
${GMQE_SCORER_DIR}/feature_groups.json
${GMQE_SCORER_DIR}/nn_0.dat
${GMQE_SCORER_DIR}/nn_1.dat
${GMQE_SCORER_DIR}/nn_2.dat
${GMQE_SCORER_DIR}/nn_3.dat
)

add_custom_target(qmean_data_files ALL)

copy_if_different("./" "${SHARED_DATA_PATH}/potentials" 
                  "${POTENTIAL_FILES}" "QMEAN_POTENTIAL_FILES" qmean_data_files)

copy_if_different("./" "${SHARED_DATA_PATH}/reference_values" 
                  "${REFERENCE_TABLE_FILES}" "QMEAN_REFERENCE_TABLE_FILES" 
                  qmean_data_files)         

copy_if_different("./" "${SHARED_DATA_PATH}/scorer" 
                  "${LINEAR_SCORER_FILES}" "QMEAN_LINEAR_SCORER_FILES" 
                  qmean_data_files)  

copy_if_different("./" "${SHARED_DATA_PATH}/scorer/local_nn_scorer" 
                  "${LOCAL_NN_SCORER_FILES}" "QMEAN_LOCAL_NN_SCORER_FILES" 
                  qmean_data_files) 

copy_if_different("./" "${SHARED_DATA_PATH}/scorer/gmqe_scorer" 
                  "${GMQE_SCORER_FILES}" "QMEAN_GMQE_SCORER_FILES" 
                  qmean_data_files)  


install(FILES ${POTENTIAL_FILES} DESTINATION "share/qmean/potentials")
install(FILES ${REFERENCE_TABLE_FILES} DESTINATION "share/qmean/reference_values")
install(FILES ${LINEAR_SCORER_FILES} DESTINATION "share/qmean/scorer")
install(FILES ${LOCAL_NN_SCORER_FILES} DESTINATION "share/qmean/scorer/local_nn_scorer")
install(FILES ${GMQE_SCORER_FILES} DESTINATION "share/qmean/scorer/gmqe_scorer")

