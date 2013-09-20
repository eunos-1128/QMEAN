#include <boost/python.hpp>
#include <boost/python/register_ptr_to_python.hpp>

#include <qmean/cbeta_statistics.hh>
#include <qmean/cbeta_potential.hh>
#include <qmean/vec_list_magic.hh>

using namespace qmean;
using namespace boost::python;

Real GetCountatomdist(CBetaStatisticPtr s, ost::conop::AminoAcid a, ost::conop::AminoAcid b, int dist_bin) { return s->GetCount(a,b,dist_bin); }
Real GetCountatom(CBetaStatisticPtr s, ost::conop::AminoAcid a, ost::conop::AminoAcid b) { return s->GetCount(a,b); }
Real GetCountdist(CBetaStatisticPtr s, int dist_bin) { return s->GetCount(dist_bin); }

Real GetEnergyatomdist(CBetaPotentialPtr p, ost::conop::AminoAcid a, ost::conop::AminoAcid b, Real dist) { return p->GetEnergy(a, b, dist); }
Real GetEnergyresidueenv(CBetaPotentialPtr p, ost::mol::ResidueView& target, ost::mol::EntityView& env, bool normalize) { return p->GetEnergy(target, env, normalize); }
Real GetEnergyres(CBetaPotentialPtr p, ost::mol::ResidueView& target, bool normalize) { return p->GetEnergy(target, normalize); }

list WrapGetEnergies(CBetaPotentialPtr p, ost::mol::EntityView& target, ost::mol::EntityView& env, bool normalize){
  std::vector<Real> energies = p->GetEnergies(target, env, normalize);
  list ret_list = VecToList<Real>(energies);
  return ret_list;
}

void export_CBeta(){

  class_<impl::CBetaOpts>("CBetaOpts", no_init)
    .def_readonly("lower_cutoff", &impl::CBetaOpts::lower_cutoff)
    .def_readonly("upper_cutoff", &impl::CBetaOpts::upper_cutoff)
    .def_readonly("number_of_bins", &impl::CBetaOpts::number_of_bins)
    .def_readonly("sequence_sep", &impl::CBetaOpts::sequence_sep)
    .def_readonly("sigma", &impl::CBetaOpts::sigma)
  ;

  class_<CBetaStatistic, bases<StatisticBase> > ("CBetaStatistic", init<Real, Real, int, int>())
    .def("Load", &CBetaStatistic::Load, (arg("filename"))).staticmethod("Load")
    .def("Save", &CBetaStatistic::Save, (arg("filename")))
    .def("Extract", &CBetaStatistic::Extract, (arg("target"),arg("env"),arg("weight")=1.0))
    .def("GetTotalCount", &CBetaStatistic::GetTotalCount)
    .def("GetCount", &GetCountatomdist, (arg("amino_acid_a"), arg("amino_acid_b"), arg("dist_bin")))
    .def("GetCount", &GetCountatom, (arg("amino_acid_a"), arg("amino_acid_b")))
    .def("GetCount", &GetCountdist, (arg("dist_bin")))
    .def("GetOptions", &CBetaStatistic::GetOpts, return_value_policy<reference_existing_object>())
  ;

  class_<CBetaPotential, bases<PotentialBase> >("CBetaPotential", no_init)
    .def("Load", &CBetaPotential::Load, (arg("filename"))).staticmethod("Load")
    .def("Save", &CBetaPotential::Save, (arg("filename")))
    .def("Create", &CBetaPotential::Create, (arg("cbeta_statistic"), arg("sigma")=0.02, arg("reference_state")="classic")).staticmethod("Create")
    .def("SetEnvironment", &CBetaPotential::SetEnvironment, args("environment"))
    .def("GetEnergy", &GetEnergyatomdist, (arg("atom_a"),arg("atom_b"), arg("dist")))
    .def("GetEnergy", &GetEnergyresidueenv, (arg("target_residue"), arg("env"), arg("normalize")=true))
    .def("GetEnergy", &GetEnergyres, (arg("target_residue"),arg("normalize")=true))
    .def("GetEnergies", &WrapGetEnergies, (arg("target"), arg("environment"),arg("normalize")=true))
    .def("GetTotalEnergy", &CBetaPotential::GetTotalEnergy, (arg("target_view"), arg("environment_view"), arg("normalize")=true))
    .def("GetCounts", &CBetaPotential::GetCounts)
    .def("GetOptions", &CBetaPotential::GetOpts, return_value_policy<reference_existing_object>())
  ;

  register_ptr_to_python<CBetaStatisticPtr>(); 
  register_ptr_to_python<CBetaPotentialPtr>();
  
}