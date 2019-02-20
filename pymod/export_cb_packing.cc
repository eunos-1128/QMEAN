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

#include <qmean/cb_packing_statistics.hh>
#include <qmean/cb_packing_potential.hh>
#include <qmean/vec_list_magic.hh>

using namespace qmean;
using namespace boost::python;

Real GetCountamino_acid(CBPackingStatisticPtr s, ost::conop::AminoAcid aa) { return s->GetCount(aa); }
Real GetCountamino_acidcount(CBPackingStatisticPtr s, ost::conop::AminoAcid aa, int count) { return s->GetCount(aa, count); }
Real GetCountcount(CBPackingStatisticPtr s, int count) { return s->GetCount(count); }

Real GetEnergy(CBPackingPotentialPtr p, ost::conop::AminoAcid aa, int count) { return p->GetEnergy(aa,count); }
Real GetEnergyResEnv(CBPackingPotentialPtr p, ost::mol::ResidueView& target, ost::mol::EntityView& env) { return p->GetEnergy(target, env); }
Real GetEnergyRes(CBPackingPotentialPtr p, ost::mol::ResidueView& target) { return p->GetEnergy(target); }

list WrapGetEnergies(CBPackingPotentialPtr p, ost::mol::EntityView& target, ost::mol::EntityView& env){
  std::vector<Real> energies = p->GetEnergies(target, env);
  list ret_list = VecToList<Real>(energies);
  return ret_list;
}

void export_CBPacking()
{
  class_<impl::CBPackingOpts>("CBPackingOpts", no_init)
    .def_readonly("cutoff", &impl::CBPackingOpts::cutoff)
    .def_readonly("max_counts", &impl::CBPackingOpts::max_counts)
    .def_readonly("sigma", &impl::CBPackingOpts::sigma)
  ;

  class_<CBPackingStatistic, bases<StatisticBase> >("CBPackingStatistic", init<Real, int>())
    .def("Load", &CBPackingStatistic::Load,(arg("filename"))).staticmethod("Load")
    .def("Save", &CBPackingStatistic::Save,(arg("filename")))
    .def("Extract", &CBPackingStatistic::Extract, (arg("target"), arg("env"), arg("weight")=1.0))  
    .def("GetTotalCount", &CBPackingStatistic::GetTotalCount)
    .def("GetCount", &GetCountcount, (arg("count")))
    .def("GetCount", &GetCountamino_acid, (arg("amino_acid")))
    .def("GetCount", &GetCountamino_acidcount, (arg("amino_acid"),arg("count")))
    .def("GetOptions", &CBPackingStatistic::GetOpts)
  ;

  
  class_<CBPackingPotential, bases<PotentialBase> >("CBPackingPotential", no_init)
    .def("Load", CBPackingPotential::Load, (arg("filename"))).staticmethod("Load")
    .def("Save", &CBPackingPotential::Save, (arg("filename")))
    .def("Create",&CBPackingPotential::Create, (arg("packing_statistic"), arg("sigma")=0.02, arg("reference_state")="classic")).staticmethod("Create")
    .def("SetEnvironment", &CBPackingPotential::SetEnvironment, (arg("environment")))
    .def("GetEnergy", &GetEnergy, (arg("amino_acid"),arg("count")))   
    .def("GetEnergy", &GetEnergyResEnv, (arg("target_residue"),arg("environment")))
    .def("GetEnergy", &GetEnergyRes, (arg("target_residue")))
    .def("GetEnergies", &WrapGetEnergies, (arg("target_view"), arg("environment_view")))
    .def("GetTotalEnergy", &CBPackingPotential::GetTotalEnergy, (arg("target_view"), arg("environment_view"),arg("normalize")=true))
    .def("GetOptions", &CBPackingPotential::GetOpts)
  ;


  register_ptr_to_python<CBPackingStatisticPtr>();
  register_ptr_to_python<CBPackingPotentialPtr>();

}

