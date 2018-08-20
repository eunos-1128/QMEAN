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


void export_Packing()
{
  class_<impl::PackingOpts>("PackingOpts", no_init)
    .def_readonly("cutoff", &impl::PackingOpts::cutoff)
    .def_readonly("max_counts", &impl::PackingOpts::max_counts)
    .def_readonly("sigma", &impl::PackingOpts::sigma)
  ;

  class_<PackingStatistic, bases<StatisticBase> >("PackingStatistic", init<Real, int>())
    .def("Load", &PackingStatistic::Load,(arg("filename"))).staticmethod("Load")
    .def("Save", &PackingStatistic::Save,(arg("filename")))
    .def("Extract", &PackingStatistic::Extract, (arg("target"), arg("env"), arg("weight")=1.0))  
    .def("GetTotalCount", &PackingStatistic::GetTotalCount)
    .def("GetCount", &GetCountcount, (arg("count")))
    .def("GetCount", &GetCountatom, (arg("atom")))
    .def("GetCount", &GetCountatomcount, (arg("atom"),arg("count")))
    .def("GetOptions", &PackingStatistic::GetOpts)
  ;

  class_<PackingPotential, bases<PotentialBase> >("PackingPotential", no_init)
    .def("Load", PackingPotential::Load, (arg("filename"))).staticmethod("Load")
    .def("Save", &PackingPotential::Save, (arg("filename")))
    .def("Create", &WrapCreate, (arg("packing_statistic"), arg("sigma")=0.02, arg("reference_state")="classic")).staticmethod("Create")    
    .def("SetEnvironment", &PackingPotential::SetEnvironment, (arg("environment")))
    .def("GetEnergy", &GetEnergy, (arg("atom"),arg("count")))   
    .def("GetEnergy", &GetEnergyResEnv, (arg("target_residue"),arg("environment"),arg("normalize")=true))
    .def("GetEnergy", &GetEnergyRes, (arg("target_residue"),arg("normalize")=true))
    .def("GetEnergies", &WrapGetEnergies, (arg("target_view"), arg("environment_view"),arg("normalize")=true))
    .def("GetTotalEnergy", &PackingPotential::GetTotalEnergy, (arg("target_view"), arg("environment_view"),arg("normalize")=true))
    .def("GetOptions", &PackingPotential::GetOpts)
  ;


  register_ptr_to_python<PackingStatisticPtr>();
  register_ptr_to_python<PackingPotentialPtr>();

}

