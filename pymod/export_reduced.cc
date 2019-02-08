// Copyright (c) 2013-2018, SIB - Swiss Institute of Bioinformatics and
// Biozentrum - University of Basel
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <boost/python.hpp>
#include <boost/python/register_ptr_to_python.hpp>

#include <qmean/reduced_potential.hh>
#include <qmean/reduced_statistics.hh>
#include <qmean/vec_list_magic.hh>


using namespace qmean;
using namespace boost::python;

Real GetCountaa(ReducedStatisticPtr s, ost::conop::AminoAcid aa_one, ost::conop::AminoAcid aa_two) { return s->GetCount(aa_one, aa_two); }
Real GetCountaabin(ReducedStatisticPtr s, ost::conop::AminoAcid aa_one,ost::conop::AminoAcid aa_two, int a, int b, int c, int d) { return s->GetCount(aa_one, aa_two, a, b, c, d); }
Real GetCountbin(ReducedStatisticPtr s, int a, int b, int c, int d){ return s->GetCount(a,b,c,d); }

Real GetEnergy(ReducedPotentialPtr p, ost::conop::AminoAcid aa_one, ost::conop::AminoAcid aa_two, Real dist, Real alpha, Real beta, Real gamma) { return p->GetEnergy(aa_one, aa_two, dist, alpha, beta, gamma); }
Real GetEnergyRes(ReducedPotentialPtr p, ost::mol::ResidueView& res, bool normalize) { p->GetEnergy(res, normalize); }
Real GetEnergyResEnv(ReducedPotentialPtr p, ost::mol::ResidueView& res, ost::mol::EntityView& env, bool normalize) { p->GetEnergy(res, env, normalize); }

list WrapGetEnergies(ReducedPotentialPtr p, ost::mol::EntityView& target, ost::mol::EntityView& env, bool normalize){
  std::vector<Real> energies = p->GetEnergies(target,env, normalize);
  list ret_list = VecToList<Real>(energies);
  return ret_list;
}

void export_Reduced()
{

  class_<impl::ReducedOpts>("ReducedOpts", no_init)
    .def_readonly("lower_cutoff", &impl::ReducedOpts::lower_cutoff)
    .def_readonly("upper_cutoff", &impl::ReducedOpts::upper_cutoff)
    .def_readonly("num_angle_bins", &impl::ReducedOpts::num_angle_bins)
    .def_readonly("num_dihedral_bins", &impl::ReducedOpts::num_dihedral_bins)
    .def_readonly("num_dist_bins", &impl::ReducedOpts::num_dist_bins)
    .def_readonly("sequence_sep", &impl::ReducedOpts::sequence_sep)
    .def_readonly("sigma", &impl::ReducedOpts::sigma)
  ;

  class_<ReducedStatistic, bases<StatisticBase> >("ReducedStatistic", init<Real, Real, int, int, int, int>())
    .def("Load", &ReducedStatistic::Load, (arg("filename"))).staticmethod("Load")
    .def("Save", &ReducedStatistic::Save, (arg("filename")))
    .def("Extract", &ReducedStatistic::Extract, (arg("target"), arg("env"), arg("weight")=1.0))
    .def("GetTotalCount", &ReducedStatistic::GetTotalCount)
    .def("GetCount", &GetCountaabin, (arg("aa_one"),arg("aa_two"),arg("dist_bin"),arg("alpha_bin"),arg("beta_bin"),arg("gamma_bin")))
    .def("GetCount", &GetCountbin, (arg("dist_bin"),arg("alpha_bin"),arg("beta_bin"),arg("gamma_bin")))
    .def("GetCount", &GetCountaa, (arg("aa_one"),arg("aa_two")))
    .def("GetOptions", &ReducedStatistic::GetOpts)
  ;

  class_<ReducedPotential, bases<PotentialBase> >("ReducedPotential", no_init)
    .def("Load", &ReducedPotential::Load, (arg("filename"))).staticmethod("Load")
    .def("Save", &ReducedPotential::Save, (arg("filename")))
    .def("Create", &ReducedPotential::Create, (arg("reduced_statistic"),arg("sigma")=0.02,arg("reference_state")="classic")).staticmethod("Create")
    .def("SetEnvironment", &ReducedPotential::SetEnvironment, args("environment"))
    .def("GetEnergy", &GetEnergy, (arg("amino_acid_a"),arg("amino_acid_b"), arg("distance"),arg("alpha"),arg("beta"),arg("gamma")))
    .def("GetEnergy", &GetEnergyRes, (arg("target_residue"),arg("normalize")=true))
    .def("GetEnergy", &GetEnergyResEnv, (arg("target_residue"),arg("environment_view"), arg("normalize")=true))
    .def("GetEnergies", &WrapGetEnergies, (arg("target_view"),arg("environment_view"),arg("normalize")=true))
    .def("GetTotalEnergy", &ReducedPotential::GetTotalEnergy, (arg("target_view"), arg("environment_view"), arg("normalize")=true))
    .def("GetCounts", &ReducedPotential::GetCounts)
    .def("GetOptions", &ReducedPotential::GetOpts)
  ; 

  register_ptr_to_python<ReducedStatisticPtr>();
  register_ptr_to_python<ReducedPotentialPtr>();
}

