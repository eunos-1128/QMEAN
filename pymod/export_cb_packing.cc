// QMEAN:
// ----------
// 
// Copyright (C) 2018 University of Basel and the Authors.
// Authors: Gabriel Studer, Pascal Benkert, Marco Biasini
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

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

