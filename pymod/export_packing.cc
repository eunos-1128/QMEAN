#include <boost/python.hpp>
#include <boost/python/register_ptr_to_python.hpp>

#include <qmean/packing_statistics.hh>
#include <qmean/packing_potential.hh>
#include <qmean/vec_list_magic.hh>


using namespace qmean;
using namespace boost::python;

Real GetCountatom(PackingStatisticPtr s, qmean::atom::ChemType a) { return s->GetCount(a); }
Real GetCountatomcount(PackingStatisticPtr s, qmean::atom::ChemType a, int count) { return s->GetCount(a, count); }
Real GetCountcount(PackingStatisticPtr s, int count) { return s->GetCount(count); }

Real GetEnergy(PackingPotentialPtr p, atom::ChemType a, int count) { return p->GetEnergy(a,count); }
Real GetEnergyResEnv(PackingPotentialPtr p, ost::mol::ResidueView& target, ost::mol::EntityView& env, bool normalize) { return p->GetEnergy(target, env, normalize); }
Real GetEnergyRes(PackingPotentialPtr p, ost::mol::ResidueView& target, bool normalize) { return p->GetEnergy(target, normalize); }

list WrapGetEnergies(PackingPotentialPtr p, ost::mol::EntityView& target, ost::mol::EntityView& env, bool normalize){
  std::vector<Real> energies = p->GetEnergies(target, env, normalize);
  list ret_list = VecToList<Real>(energies);
  return ret_list;
}

PackingPotentialPtr WrapCreate(PackingStatisticPtr s, Real sigma, const String& reference){
  PackingPotentialPtr p = PackingPotential::Create(s,sigma,reference);
  return p;
}

PackingPotentialPtr WrapHybridCreate(PackingStatisticPtr s1, PackingStatisticPtr s2, Real sigma, const String& reference, Real max_energy){
  PackingPotentialPtr p = PackingPotential::Create(s1,s2,sigma,reference,max_energy);
  return p;
}

void export_Packing()
{
  class_<impl::PackingOpts>("PackingOpts", no_init)
    .def_readonly("cutoff", &impl::PackingOpts::cutoff)
    .def_readonly("max_counts", &impl::PackingOpts::max_counts)
    .def_readonly("bin_size", &impl::PackingOpts::bin_size)
    .def_readonly("sigma", &impl::PackingOpts::sigma)
  ;

  class_<PackingStatistic, bases<StatisticBase> >("PackingStatistic", init<Real, int, int>())
    .def("Load", &PackingStatistic::Load,(arg("filename"))).staticmethod("Load")
    .def("Save", &PackingStatistic::Save,(arg("filename")))
    .def("Extract", &PackingStatistic::Extract, (arg("target"), arg("env"), arg("weight")=1.0))  
    .def("GetTotalCount", &PackingStatistic::GetTotalCount)
    .def("GetCount", &GetCountcount, (arg("count")))
    .def("GetCount", &GetCountatom, (arg("atom")))
    .def("GetCount", &GetCountatomcount, (arg("atom"),arg("count")))
    .def("GetOptions", &PackingStatistic::GetOpts, return_value_policy<reference_existing_object>())
  ;

  class_<PackingPotential, bases<PotentialBase> >("PackingPotential", no_init)
    .def("Load", PackingPotential::Load, (arg("filename"))).staticmethod("Load")
    .def("Save", &PackingPotential::Save, (arg("filename")))
    .def("Create", &WrapHybridCreate, (arg("packing_statistic_one"),arg("packing_statistic_two"),arg("sigma")=0.02,arg("reference_state")="classic",arg("max_energy")=8.0))
    .def("Create", &WrapCreate, (arg("packing_statistic"), arg("sigma")=0.02, arg("reference_state")="classic")).staticmethod("Create")    
    .def("SetEnvironment", &PackingPotential::SetEnvironment, (arg("environment")))
    .def("GetEnergy", &GetEnergy, (arg("atom"),arg("count")))   
    .def("GetEnergy", &GetEnergyResEnv, (arg("target_residue"),arg("environment"),arg("normalize")=true))
    .def("GetEnergy", &GetEnergyRes, (arg("target_residue"),arg("normalize")=true))
    .def("GetEnergies", &WrapGetEnergies, (arg("target_view"), arg("environment_view"),arg("normalize")=true))
    .def("GetTotalEnergy", &PackingPotential::GetTotalEnergy, (arg("target_view"), arg("environment_view"),arg("normalize")=true))
    .def("GetOptions", &PackingPotential::GetOpts, return_value_policy<reference_existing_object>())
  ;


  register_ptr_to_python<PackingStatisticPtr>();
  register_ptr_to_python<PackingPotentialPtr>();

}

