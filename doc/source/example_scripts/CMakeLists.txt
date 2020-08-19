set(DOC_UNIT_TESTS
  test_doctests.py
)

set(DOC_TEST_DATA
  example_data/1CRN.pdb
  example_data/7ODC.pdb
  example_data/2QUD.pdb 
  example_data/3T47.pdb 
  example_data/4C1A.pdb 
  example_data/3ZSJ.pdb 
  example_data/4I4T.pdb 
  example_data/3H6R.pdb 
  example_data/3BSO.pdb                
  example_data/3M0Z.pdb 
  example_data/2GPI.pdb 
  example_data/3MHP.pdb 
  example_data/3V19.pdb 
  example_data/2IXS.pdb 
  example_data/4G6T.pdb 
  example_data/3ZG9.pdb 
  example_data/3WSG.pdb                  
  example_data/1JH6.pdb 
  example_data/2XZE.pdb 
  example_data/2OLR.pdb 
  example_data/1TIF.pdb 
  example_data/3C4S.pdb 
  example_data/2CVI.pdb 
  example_data/4EET.pdb 
  example_data/1LM5.pdb                  
  example_data/1MJ5.pdb 
  example_data/2FJR.pdb 
  example_data/2D3G.pdb 
  example_data/3ZNV.pdb 
  example_data/3WA2.pdb 
  example_data/3WU2.pdb 
  example_data/4M7R.pdb 
  example_data/2PR7.pdb                  
  example_data/3FTJ.pdb 
  example_data/1KD8.pdb 
  example_data/1HBN.pdb 
  example_data/4TRK.pdb 
  example_data/2GB4.pdb 
  example_data/3HNY.pdb 
  example_data/4R7Q.pdb 
  example_data/1EAQ.pdb                  
  example_data/2O1Q.pdb 
  example_data/4DX5.pdb 
  example_data/1XAK.pdb 
  example_data/5CSM.pdb 
  example_data/2XWP.pdb 
  example_data/2UWA.pdb 
  example_data/3SQZ.pdb 
  example_data/4MT8.pdb                  
  example_data/3HTU.pdb
  example_data/3EML.pdb
  example_data/1pvv.1.A.fasta 
  example_data/1v1v.1.A.fasta 
  example_data/2w37.1.A.fasta 
  example_data/3grf.1.A.fasta
  example_data/4nf2.1.A.fasta
  example_data/1pvv.1.A.pdb 
  example_data/1v1v.1.A.pdb
  example_data/2w37.1.A.pdb 
  example_data/3grf.1.A.pdb
  example_data/4nf2.1.A.pdb
  example_data/target_seq.fasta
  example_data/model_01.pdb
  example_data/original_hhblits_alignment.pdb
  example_data/shift_in_front_helix_four.pdb
  example_data/shift_into_middle.pdb
  example_data/shift_towards_cter.pdb
  example_data/original_hhblits_alignment_local_scores.txt
  example_data/shift_in_front_helix_four_local_scores.txt
  example_data/shift_into_middle_local_scores.txt
  example_data/shift_towards_cter_local_scores.txt
  example_data/housing.csv
  example_data/1crn.1.A.hhm
  example_data/1crn_dc.dat
  example_data/1crn_3szs_aln.fasta
  example_data/3szs.2.A.pdb
  example_data/3szs.2.A.hhm
)

set (DOC_TEST_SCRIPTS
  torsion_potential_example.py
  ss_agreement_example.py
  smoother_example.py
  psipred_handler_example.py
  packing_potential_example.py
  membrane_example.py
  local_scorer_example.py
  interaction_potential_example.py
  disco_example.py
  container_example.py
  chemtype_example.py
  cbeta_potential_example.py
  accpro_handler_example.py
  assess_model_quality_example.py
  assess_membrane_model_quality.py
  reproduce_fig4_from_publication.py
  regressor_training.py
  gmqe_example.py
)

add_custom_target(doctest)
add_dependencies(check doctest)
add_dependencies(doctest qmean)

qmean_unittest(MODULE doc
  SOURCES "${DOC_UNIT_TESTS}"
  DATA "${DOC_TEST_DATA};${DOC_TEST_SCRIPTS}"
  BASE_TARGET doctest
)
