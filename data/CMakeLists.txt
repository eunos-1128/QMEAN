set(POTENTIAL_FILES
qmean/potentials/soluble_global_potentials.dat 
qmean/potentials/soluble_local_potentials.dat 
qmean/potentials/membrane_local_potentials.dat
qmean/potentials/membrane_global_potentials.dat
)
set(REFERENCES
qmean/reference_values/reference_tab_promod_scorer.txt
)
set(SCORER
qmean/scorer/promod_local_scorer.dat
qmean/scorer/promod_global_scorer.dat
qmean/scorer/membrane_local_scorer.dat
qmean/scorer/disco_tree.pkl
)

add_custom_target(qmean_data_files ALL)
copy_if_different("./" "${SHARED_DATA_PATH}/potentials" 
                  "${POTENTIAL_FILES}" "qmean_DATA_FILES" qmean_data_files)
copy_if_different("./" "${SHARED_DATA_PATH}/references" 
                  "${REFERENCES}" "REFERENCES" qmean_data_files)         
copy_if_different("./" "${SHARED_DATA_PATH}/scorer" 
                  "${SCORER}" "SCORER" qmean_data_files)  
install(FILES ${POTENTIAL_FILES} DESTINATION "share/qmean/potentials")
install(FILES ${REFERENCES} DESTINATION "share/qmean/references")
install(FILES ${SCORER} DESTINATION "share/qmean/scorer")
