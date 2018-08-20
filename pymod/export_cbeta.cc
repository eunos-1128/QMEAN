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

CBetaPotentialPtr WrapCreate(CBetaStatisticPtr s, Real sigma, const String& reference){
  CBetaPotentialPtr p = CBetaPotential::Create(s,sigma,reference);
  return p;
}

CBetaPotentialPtr WrapHybridCreate(CBetaStatisticPtr s1, CBetaStatisticPtr s2, Real sigma, const String& reference, Real max_energy){
  CBetaPotentialPtr p = CBetaPotential::Create(s1,s2,sigma,reference,max_energy);
  return p;
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
    .def("GetOptions", &CBetaStatistic::GetOpts)
  ;

  class_<CBetaPotential, bases<PotentialBase> >("CBetaPotential", no_init)
    .def("Load", &CBetaPotential::Load, (arg("filename"))).staticmethod("Load")
    .def("Save", &CBetaPotential::Save, (arg("filename")))
    .def("Create", &WrapHybridCreate, (arg("cbeta_statistic_one"),arg("cbeta_statistic_two"),arg("sigma")=0.02,arg("reference_state")="classic",arg("max_energy")=8.0))
    .def("Create", &WrapCreate, (arg("cbeta_statistic"), arg("sigma")=0.02, arg("reference_state")="classic")).staticmethod("Create")
    .def("SetEnvironment", &CBetaPotential::SetEnvironment, args("environment"))
    .def("GetEnergy", &GetEnergyatomdist, (arg("atom_a"),arg("atom_b"), arg("dist")))
    .def("GetEnergy", &GetEnergyresidueenv, (arg("target_residue"), arg("env"), arg("normalize")=true))
    .def("GetEnergy", &GetEnergyres, (arg("target_residue"),arg("normalize")=true))
    .def("GetEnergies", &WrapGetEnergies, (arg("target"), arg("environment"),arg("normalize")=true))
    .def("GetTotalEnergy", &CBetaPotential::GetTotalEnergy, (arg("target_view"), arg("environment_view"), arg("normalize")=true))
    .def("GetCounts", &CBetaPotential::GetCounts)
    .def("GetOptions", &CBetaPotential::GetOpts)
  ;

  register_ptr_to_python<CBetaStatisticPtr>(); 
  register_ptr_to_python<CBetaPotentialPtr>();
  
}

