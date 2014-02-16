set(QMEAN_MODULES
  __init__.py
  predicted_sequence_features.py
  conf.py
  pisces.py
  mqa.py
  mqa_result.py
  mqa_membrane.py
  mqa_result_membrane.py
  crap_checker.py
  score_calculator.py
  reference_set.py

)
set(QMEAN_PYMOD
export_interaction.cc
export_torsion.cc
export_packing.cc
export_reduced.cc
export_cbeta.cc
export_ss_agreement.cc
export_base.cc
export_spherical_smoother.cc
export_membrane.cc
wrap_qmean.cc
)

pymod(NAME qmean OUTPUT_DIR qmean CPP ${QMEAN_PYMOD} PY ${QMEAN_MODULES}
      PREFIX "")
