#include <boost/python.hpp>
#include <boost/python/register_ptr_to_python.hpp>

#include <qmean/disco.hh>
#include <qmean/vec_list_magic.hh>

using namespace qmean;
using namespace boost::python;

namespace{

void WrapAddData(DisCoContainerPtr container, 
                 const ost::seq::AlignmentHandle& aln, 
                 const geom::Vec3List& positions,
                 const list& pos_seqres_mapping) {
  std::vector<uint> v_pos_seqres_mapping = 
  ListToVec<uint>(pos_seqres_mapping);

  container->AddData(aln, positions, v_pos_seqres_mapping);
}


list WrapGetConstraint(DisCoContainerPtr container,
                       uint rnum_one, uint rnum_two) {

  const std::vector<Real>& constraint = 
  container->GetConstraint(rnum_one, rnum_two);

  return VecToList<Real>(constraint);
}


list WrapGetScores(DisCoContainerPtr container,
                   ost::mol::EntityView& view) {
  std::vector<Real> scores = container->GetScores(view);
  return VecToList<Real>(scores);
}

} // ns



void export_disco() {

  class_<DisCoContainer, DisCoContainerPtr> ("DisCoContainer", init<const ost::seq::SequenceHandle&>())
    .def("Load", &DisCoContainer::Load, (arg("filename"))).staticmethod("Load")
    .def("Save", &DisCoContainer::Save, (arg("filename")))
    .def("AddData", &WrapAddData, (arg("aln"), arg("positions"), 
                                   arg("pos_seqres_mapping")))
    .def("CalculateConstraints", &DisCoContainer::CalculateConstraints,
         (arg("dist_cutoff")=15.0, arg("gamma")=70.0, 
          arg("seqsim_clustering_cutoff")=0.4, arg("bin_size")=1.0))
    .def("HasConstraint", &DisCoContainer::HasConstraint, (arg("rnum_one"),
                                                           arg("rnum_two")))
    .def("GetConstraint", &WrapGetConstraint, (arg("rnum_one"),
                                               arg("rnum_two")))
    .def("GetScores", &WrapGetScores, (arg("view")))
    .def("GetDistCutoff", &DisCoContainer::GetDistCutoff)
    .def("GetGamma", &DisCoContainer::GetGamma)
    .def("GetSeqSimClusteringCutoff", &DisCoContainer::GetSeqSimClusteringCutoff)
    .def("GetNumBins", &DisCoContainer::GetNumBins)
    .def("GetBinSize", &DisCoContainer::GetBinSize)
  ;

}
