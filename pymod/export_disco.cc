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
  std::vector<Real> scores;
  container->FillData(view, scores);
  return VecToList<Real>(scores);
}

dict WrapGetScoresWithData(DisCoContainerPtr container,
                           ost::mol::EntityView& view) {
  std::vector<Real> scores;
  std::vector<uint> counts;
  std::vector<Real> avg_num_clusters;
  std::vector<Real> avg_max_seqsim;
  std::vector<Real> avg_max_seqid;
  std::vector<Real> avg_variance;
  container->FillData(view, scores, counts, avg_num_clusters, avg_max_seqsim, 
                      avg_max_seqid, avg_variance);

  dict return_dict;
  return_dict["scores"] = VecToList<Real>(scores);
  return_dict["counts"] = VecToList<uint>(counts);
  return_dict["avg_num_clusters"] = VecToList<Real>(avg_num_clusters);
  return_dict["avg_max_seqsim"] = VecToList<Real>(avg_max_seqsim);
  return_dict["avg_max_seqid"] = VecToList<Real>(avg_max_seqid);
  return_dict["avg_variance"] = VecToList<Real>(avg_variance);  
  return return_dict;
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
          arg("seqsim_clustering_cutoff")=0.5, arg("bin_size")=1.0))
    .def("HasConstraint", &DisCoContainer::HasConstraint, (arg("rnum_one"),
                                                           arg("rnum_two")))
    .def("GetConstraint", &WrapGetConstraint, (arg("rnum_one"),
                                               arg("rnum_two")))
    .def("GetScores", &WrapGetScores, (arg("view")))
    .def("GetScoresWithData", &WrapGetScoresWithData, (arg("view")))
    .def("GetSeqres", &DisCoContainer::GetSeqres)
    .def("GetDistCutoff", &DisCoContainer::GetDistCutoff)
    .def("GetGamma", &DisCoContainer::GetGamma)
    .def("GetSeqSimClusteringCutoff", &DisCoContainer::GetSeqSimClusteringCutoff)
    .def("GetNumBins", &DisCoContainer::GetNumBins)
    .def("GetBinSize", &DisCoContainer::GetBinSize)
  ;

}
