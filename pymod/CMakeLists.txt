set(QMEAN_MODULES
  __init__.py
  predicted_sequence_features.py
  conf.py
  mqa.py
  mqa_result.py
  mqa_membrane.py
  mqa_result_membrane.py
  score_calculator.py
  reference_set.py
  mlp_regressor.py
  gmqe.py
)

set(QMEAN_PYMOD
export_interaction.cc
export_torsion.cc
export_packing.cc
export_reduced.cc
export_cbeta.cc
export_cb_packing.cc
export_hbond.cc
export_ss_agreement.cc
export_base.cc
export_spherical_smoother.cc
export_disco.cc
export_clash.cc
export_gmqe.cc
export_extract_data_helper.cc
wrap_qmean.cc
)

pymod(NAME qmean OUTPUT_DIR qmean CPP ${QMEAN_PYMOD} PY ${QMEAN_MODULES}
      PREFIX "")

