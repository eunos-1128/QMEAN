#include <boost/python.hpp>
#include <boost/python/register_ptr_to_python.hpp>

#include <qmean/torsion_potential.hh>
#include <qmean/torsion_statistics.hh>
#include <qmean/vec_list_magic.hh>

using namespace qmean;
using namespace boost::python;

TorsionStatisticPtr WrapTorsionStatisticConstructor(boost::python::list g_id, boost::python::list a_b){
  std::vector<String> group_ids = ListToVec<String>(g_id);
  std::vector<int> angle_bins = ListToVec<int>(a_b);
  TorsionStatisticPtr s(new TorsionStatistic(group_ids, angle_bins));
  return s;
}

Real GetTotalCountWithout_gi(TorsionStatisticPtr s) { return s->GetTotalCount(); }
Real GetTotalCountWith_gi(TorsionStatisticPtr s, const String& gi){ return s->GetTotalCount(gi); }

Real GetCount_gi_bins(TorsionStatisticPtr s, const String& gi, boost::python::list& bins) {
  std::vector<int> vec_bins = ListToVec<int>(bins); 
  return s->GetCount(gi, vec_bins); 
}

Real GetCount_bins(TorsionStatisticPtr s, boost::python::list& bins) { 
  std::vector<int> vec_bins = ListToVec<int>(bins);
  return s->GetCount(vec_bins); 
}

Real GetEnergyView(TorsionPotentialPtr p, ost::mol::ResidueView& view){ return p->GetEnergy(view); }

Real GetEnergyParams(TorsionPotentialPtr p, list& r_n, list& a){ 
  std::vector<String> residue_names = ListToVec<String>(r_n);
  std::vector<Real> angles = ListToVec<Real>(a);
  return p->GetEnergy(residue_names, angles); 
}

Real GetEnergyStringID(TorsionPotentialPtr p, const String& g_id, list& a){ 
  std::vector<Real> angles = ListToVec<Real>(a);
  return p->GetEnergy(g_id, angles); 
}

list WrapGetEnergies(TorsionPotentialPtr p, ost::mol::EntityView& target){
  std::vector<Real> energies = p->GetEnergies(target);
  list ret_list = VecToList<Real>(energies);
  return ret_list;
}

void export_Torsion()
{

class_<impl::TorsionOpts>("TorsionOpts", no_init)
  .def_readonly("sigma", &impl::TorsionOpts::sigma)
  .def_readonly("number_of_bins", &impl::TorsionOpts::num_of_bins)
  .def_readonly("group_identifier", &impl::TorsionOpts::group_identifier)
;

class_<TorsionStatistic, bases<StatisticBase> >("TorsionStatistic")
  .def("__init__", (arg("group_ids"),arg("number_of_bins"),make_constructor(&WrapTorsionStatisticConstructor)))
  .def("Load", &TorsionStatistic::Load,(arg("filename"))).staticmethod("Load")
  .def("Save", &TorsionStatistic::Save, (arg("filename")))
  .def("Extract", &TorsionStatistic::Extract, (arg("target"), arg("weight")=1.0))
  .def("GetTotalCount", &GetTotalCountWith_gi)
  .def("GetTotalCount", &GetTotalCountWithout_gi, (arg("group_identifier")))
  .def("GetCount", &GetCount_gi_bins, (arg("group_identifier"),arg("bins")))
  .def("GetCount", &GetCount_bins, (arg("bins")))
  .def("GetOptions", &TorsionStatistic::GetOpts, return_value_policy<reference_existing_object>())
;

class_<TorsionPotential, bases<PotentialBase> >("TorsionPotential", no_init)
  .def("Load", &TorsionPotential::Load, (arg("filename"))).staticmethod("Load")
  .def("Save", &TorsionPotential::Save, (arg("filename")))
  .def("Create", &TorsionPotential::Create, (arg("torsion_statistic"),arg("sigma")=0.02,arg("reference_state")="classic")).staticmethod("Create")
  .def("GetEnergies", &WrapGetEnergies, (arg("target")))
  .def("GetEnergy", &GetEnergyView, (arg("target_residue")))
  .def("GetEnergy", &GetEnergyParams, (arg("residue_names"),arg("angles")))
  .def("GetEnergy", &GetEnergyStringID, (arg("g_id"),arg("angles")))
  .def("GetTotalEnergy", &TorsionPotential::GetTotalEnergy, (arg("target_view"),arg("normalize")=true))
  .def("GetOptions", &TorsionPotential::GetOpts, return_value_policy<reference_existing_object>())
  ; 

 register_ptr_to_python<TorsionStatisticPtr>();
 register_ptr_to_python<TorsionPotentialPtr>();
}

