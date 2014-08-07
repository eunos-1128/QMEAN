#include <boost/python.hpp>
#include <boost/python/register_ptr_to_python.hpp>

#include <qmean/cb_packing_statistics.hh>
#include <qmean/cb_packing_potential.hh>
#include <qmean/vec_list_magic.hh>

using namespace qmean;
using namespace boost::python;

Real GetCountamino_acid(CBPackingStatisticPtr s, ost::conop::AminoAcid aa) { return s->GetCount(aa); }
Real GetCountamino_acidcount(CBPackingStatisticPtr s, ost::conop::AminoAcid aa, int count) { return s->GetCount(aa, count); }
Real GetCountcount(CBPackingStatisticPtr s, int count) { return s->GetCount(count); }

Real GetEnergy(CBPackingPotentialPtr p, ost::conop::AminoAcid aa, int count) { return p->GetEnergy(aa,count); }
Real GetEnergyResEnv(CBPackingPotentialPtr p, ost::mol::ResidueView& target, ost::mol::EntityView& env, bool normalize) { return p->GetEnergy(target, env, normalize); }
Real GetEnergyRes(CBPackingPotentialPtr p, ost::mol::ResidueView& target, bool normalize) { return p->GetEnergy(target, normalize); }

list WrapGetEnergies(CBPackingPotentialPtr p, ost::mol::EntityView& target, ost::mol::EntityView& env, bool normalize){
  std::vector<Real> energies = p->GetEnergies(target, env, normalize);
  list ret_list = VecToList<Real>(energies);
  return ret_list;
}

void export_CBPacking()
{
  class_<impl::CBPackingOpts>("CBPackingOpts", no_init)
    .def_readonly("cutoff", &impl::CBPackingOpts::cutoff)
    .def_readonly("max_counts", &impl::CBPackingOpts::max_counts)
    .def_readonly("bin_size", &impl::CBPackingOpts::bin_size)
    .def_readonly("sigma", &impl::CBPackingOpts::sigma)
  ;

  class_<CBPackingStatistic, bases<StatisticBase> >("CBPackingStatistic", init<Real, int, int>())
    .def("Load", &CBPackingStatistic::Load,(arg("filename"))).staticmethod("Load")
    .def("Save", &CBPackingStatistic::Save,(arg("filename")))
    .def("Extract", &CBPackingStatistic::Extract, (arg("target"), arg("env"), arg("weight")=1.0))  
    .def("GetTotalCount", &CBPackingStatistic::GetTotalCount)
    .def("GetCount", &GetCountcount, (arg("count")))
    .def("GetCount", &GetCountamino_acid, (arg("amino_acid")))
    .def("GetCount", &GetCountamino_acidcount, (arg("amino_acid"),arg("count")))
    .def("GetOptions", &CBPackingStatistic::GetOpts, return_value_policy<reference_existing_object>())
  ;

  
  class_<CBPackingPotential, bases<PotentialBase> >("CBPackingPotential", no_init)
    .def("Load", CBPackingPotential::Load, (arg("filename"))).staticmethod("Load")
    .def("Save", &CBPackingPotential::Save, (arg("filename")))
    .def("Create",&CBPackingPotential::Create, (arg("packing_statistic"), arg("sigma")=0.02, arg("reference_state")="classic")).staticmethod("Create")
    .def("SetEnvironment", &CBPackingPotential::SetEnvironment, (arg("environment")))
    .def("GetEnergy", &GetEnergy, (arg("amino_acid"),arg("count")))   
    .def("GetEnergy", &GetEnergyResEnv, (arg("target_residue"),arg("environment"),arg("normalize")=true))
    .def("GetEnergy", &GetEnergyRes, (arg("target_residue"),arg("normalize")=true))
    .def("GetEnergies", &WrapGetEnergies, (arg("target_view"), arg("environment_view"),arg("normalize")=true))
    .def("GetTotalEnergy", &CBPackingPotential::GetTotalEnergy, (arg("target_view"), arg("environment_view"),arg("normalize")=true))
    .def("GetOptions", &CBPackingPotential::GetOpts, return_value_policy<reference_existing_object>())
  ;


  register_ptr_to_python<CBPackingStatisticPtr>();
  register_ptr_to_python<CBPackingPotentialPtr>();

}

