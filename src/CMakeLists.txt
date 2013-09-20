set(QMEAN_HEADERS
module_config.hh
histogram.hh
float_histogram.hh
index.hh
multi_classifier.hh
atom_types.hh
statistic_base.hh
potential_base.hh
torsion_statistics.hh
torsion_potential.hh
interaction_statistics.hh
interaction_potential.hh
packing_statistics.hh
packing_potential.hh
reduced_statistics.hh
reduced_potential.hh
cbeta_statistics.hh
cbeta_potential.hh
ss_agreement.hh
spherical_smoother.hh
vec_list_magic.hh
membrane.hh
impl/reduced_impl.hh
impl/torsion_impl.hh
impl/packing_impl.hh
impl/interaction_impl.hh
impl/cbeta_impl.hh
)

set(QMEAN_SOURCES
atom_types.cc
statistic_base.cc
potential_base.cc
reduced_statistics.cc
reduced_potential.cc
torsion_statistics.cc
torsion_potential.cc
interaction_statistics.cc
interaction_potential.cc
packing_statistics.cc
packing_potential.cc
cbeta_statistics.cc
cbeta_potential.cc
ss_agreement.cc
spherical_smoother.cc
membrane.cc
impl/reduced_impl.cc
impl/torsion_impl.cc
impl/packing_impl.cc
impl/interaction_impl.cc
impl/cbeta_impl.cc
)


module(NAME qmean SOURCES ${QMEAN_SOURCES}
       HEADERS reduced_impl.hh IN_DIR impl ${QMEAN_HEADERS} 
       HEADER_OUTPUT_DIR qmean PREFIX ""
       LINK ${OST_LIBRARIES} ${BOOST_LIBRARIES})
