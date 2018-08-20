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

#include <qmean/hbond_statistics.hh>
#include <qmean/hbond_potential.hh>
#include <qmean/vec_list_magic.hh>

using namespace qmean;
using namespace boost::python;

Real GetCountstate(HBondStatisticPtr s, int state) { return s->GetCount(state); }
Real GetCountallinfo(HBondStatisticPtr s, int state, int d_bin, int alpha_bin,
                     int beta_bin, int gamma_bin) { 
  return s->GetCount(state, d_bin, alpha_bin, beta_bin, gamma_bin); 
}

Real GetEnergy(HBondPotentialPtr p, int state, Real d, Real alpha, Real beta, Real gamma){ 
  return p->GetEnergy(state,d,alpha,beta,gamma); 
}

Real GetEnergyResEnv(HBondPotentialPtr p, ost::mol::ResidueView& target, 
                     ost::mol::EntityView& env){ 
  return p->GetEnergy(target, env); 
}

Real GetEnergyRes(HBondPotentialPtr p, ost::mol::ResidueView& target) { 
  return p->GetEnergy(target); 
}

list WrapGetEnergies(HBondPotentialPtr p, ost::mol::EntityView& view){
  std::vector<Real> energies = p->GetEnergies(view);
  list ret_list = VecToList<Real>(energies);
  return ret_list;
}

void export_HBond()
{
  class_<impl::HBondOpts>("HBondOpts", no_init)
    .def_readonly("d_min", &impl::HBondOpts::d_min)
    .def_readonly("d_max", &impl::HBondOpts::d_max)
    .def_readonly("alpha_min", &impl::HBondOpts::alpha_min)
    .def_readonly("alpha_max", &impl::HBondOpts::alpha_max)
    .def_readonly("beta_min", &impl::HBondOpts::beta_min)
    .def_readonly("beta_max", &impl::HBondOpts::beta_max)
    .def_readonly("gamma_min", &impl::HBondOpts::gamma_min)
    .def_readonly("gamma_max", &impl::HBondOpts::gamma_max)
    .def_readonly("d_bins",&impl::HBondOpts::d_bins)
    .def_readonly("alpha_bins",&impl::HBondOpts::alpha_bins)
    .def_readonly("beta_bins",&impl::HBondOpts::beta_bins)
    .def_readonly("gamma_bins",&impl::HBondOpts::gamma_bins)
    .def_readonly("seq_sep",&impl::HBondOpts::seq_sep)
  ;

  class_<HBondStatistic, bases<StatisticBase> >("HBondStatistic", init<Real, Real, Real, Real, 
                                                                       Real, Real, Real, Real,
                                                                       int, int, int, int, int>())
    .def("Load", &HBondStatistic::Load,(arg("filename"))).staticmethod("Load")
    .def("Save", &HBondStatistic::Save,(arg("filename")))
    .def("Extract", &HBondStatistic::Extract, (arg("view"), arg("weight")=1.0))  
    .def("GetTotalCount", &HBondStatistic::GetTotalCount)
    .def("GetCount", &GetCountstate, (arg("state")))
    .def("GetCount", &GetCountallinfo, (arg("state"),arg("d_bin"),arg("alpha_bin"),
                                        arg("beta_bin"),arg("gamma_bin")))
    .def("GetOptions", &HBondStatistic::GetOpts)
  ;

  
  class_<HBondPotential, bases<PotentialBase> >("HBondPotential", no_init)
    .def("Load",&HBondPotential::Load, (arg("filename"))).staticmethod("Load")
    .def("Save", &HBondPotential::Save, (arg("filename")))
    .def("Create",&HBondPotential::Create, (arg("hbond_statistic"))).staticmethod("Create")
    .def("SetEnvironment", &HBondPotential::SetEnvironment, (arg("environment")))
    .def("GetEnergy", &GetEnergy, (arg("state"),arg("d"),arg("alpha"),arg("beta"),arg("gamma")))   
    .def("GetEnergy", &GetEnergyResEnv, (arg("target_residue"),arg("environment")))
    .def("GetEnergy", &GetEnergyRes, (arg("target_residue")))
    .def("GetEnergies", &WrapGetEnergies, (arg("target_view")))
    .def("GetTotalEnergy", &HBondPotential::GetTotalEnergy, (arg("view"),arg("normalize")=true))
    .def("GetOptions", &HBondPotential::GetOpts)
  ;

  register_ptr_to_python<HBondStatisticPtr>();
  register_ptr_to_python<HBondPotentialPtr>();
}

