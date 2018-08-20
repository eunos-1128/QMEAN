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

#include <qmean/interaction_statistics.hh>
#include <qmean/interaction_potential.hh>
#include <qmean/vec_list_magic.hh>


using namespace qmean; 
using namespace boost::python;

Real GetCountatomdist(InteractionStatisticPtr s, atom::ChemType a, atom::ChemType b, int dist_bin) { return s->GetCount(a,b,dist_bin); }
Real GetCountatom(InteractionStatisticPtr s, atom::ChemType a, atom::ChemType b) { return s->GetCount(a,b); }
Real GetCountdist(InteractionStatisticPtr s, int dist_bin) { return s->GetCount(dist_bin); }
  
Real GetEnergyatomdist(InteractionPotentialPtr p, atom::ChemType a, atom::ChemType b, Real dist) { return p->GetEnergy(a, b, dist); }
Real GetEnergyresidueenv(InteractionPotentialPtr p, ost::mol::ResidueView& target, ost::mol::EntityView& env, bool normalize) { return p->GetEnergy(target, env, normalize); }
Real GetEnergyres(InteractionPotentialPtr p, ost::mol::ResidueView& target, bool normalize) { return p->GetEnergy(target, normalize); }

list WrapGetEnergies(InteractionPotentialPtr p, ost::mol::EntityView& target, ost::mol::EntityView& env, bool normalize){
  std::vector<Real> energies = p->GetEnergies(target, env, normalize);
  list ret_list=VecToList<Real>(energies);
  return ret_list;
}

InteractionPotentialPtr WrapCreate(InteractionStatisticPtr s, Real sigma, const String& reference){
  InteractionPotentialPtr p = InteractionPotential::Create(s,sigma,reference);
  return p;
}

InteractionPotentialPtr WrapHybridCreate(InteractionStatisticPtr s1, InteractionStatisticPtr s2, Real sigma, const String& reference, Real max_energy){
  InteractionPotentialPtr p = InteractionPotential::Create(s1,s2,sigma,reference,max_energy);
  return p;
}


void export_Interaction()
{
  enum_<qmean::atom::ChemType>("ChemType")
    .value("O_GLY", atom::O_GLY)
    .value("N_GLY", atom::N_GLY)
    .value("C_GLY", atom::C_GLY)    
    .value("C_GLY_A", atom::C_GLY_A)

    .value("O_ALA", atom::O_ALA)
    .value("N_ALA", atom::N_ALA)
    .value("C_ALA", atom::C_ALA)    
    .value("C_ALA_A", atom::C_ALA_A)
    .value("C_ALA_B", atom::C_ALA_B)

    .value("O_VAL", atom::O_VAL)
    .value("N_VAL", atom::N_VAL)
    .value("C_VAL", atom::C_VAL)    
    .value("C_VAL_A", atom::C_VAL_A)
    .value("C_VAL_B", atom::C_VAL_B)
    .value("C_VAL_G", atom::C_VAL_G)

    .value("O_LEU", atom::O_LEU)
    .value("N_LEU", atom::N_LEU)
    .value("C_LEU", atom::C_LEU)  
    .value("C_LEU_A", atom::C_LEU_A)
    .value("C_LEU_B", atom::C_LEU_B)
    .value("C_LEU_G", atom::C_LEU_G)
    .value("C_LEU_D", atom::C_LEU_D)  

    .value("O_ILE", atom::O_ILE)
    .value("N_ILE", atom::N_ILE)
    .value("C_ILE", atom::C_ILE)
    .value("C_ILE_A", atom::C_ILE_A)
    .value("C_ILE_B", atom::C_ILE_B)
    .value("C_ILE_G1", atom::C_ILE_G1)
    .value("C_ILE_G2", atom::C_ILE_G2)
    .value("C_ILE_D1", atom::C_ILE_D1)

    .value("O_THR", atom::O_THR)
    .value("N_THR", atom::N_THR)
    .value("C_THR", atom::C_THR)
    .value("C_THR_A", atom::C_THR_A)
    .value("C_THR_B", atom::C_THR_B)
    .value("O_THR_G1", atom::O_THR_G1)
    .value("C_THR_G2", atom::C_THR_G2)

    .value("O_SER", atom::O_SER)
    .value("N_SER", atom::N_SER)
    .value("C_SER", atom::C_SER)
    .value("C_SER_A", atom::C_SER_A)
    .value("C_SER_B", atom::C_SER_B)
    .value("O_SER_G", atom::O_SER_G)

    .value("O_CYS", atom::O_CYS)
    .value("N_CYS", atom::N_CYS)
    .value("C_CYS", atom::C_CYS)
    .value("C_CYS_A", atom::C_CYS_A)
    .value("C_CYS_B", atom::C_CYS_B)
    .value("S_CYS_G", atom::S_CYS_G)    

    .value("O_MET", atom::O_MET)
    .value("N_MET", atom::N_MET)
    .value("C_MET", atom::C_MET)
    .value("C_MET_A", atom::C_MET_A)
    .value("C_MET_B", atom::C_MET_B)
    .value("C_MET_G", atom::C_MET_G)
    .value("S_MET_D", atom::S_MET_D)
    .value("C_MET_E", atom::C_MET_E)

    .value("O_ASP", atom::O_ASP)
    .value("N_ASP", atom::N_ASP)
    .value("C_ASP", atom::C_ASP)
    .value("C_ASP_A", atom::C_ASP_A)
    .value("C_ASP_B", atom::C_ASP_B)
    .value("C_ASP_G", atom::C_ASP_G)
    .value("O_ASP_D", atom::O_ASP_D)

    .value("O_GLU", atom::O_GLU)
    .value("N_GLU", atom::N_GLU)
    .value("C_GLU", atom::C_GLU)
    .value("C_GLU_A", atom::C_GLU_A)
    .value("C_GLU_B", atom::C_GLU_B)
    .value("C_GLU_G", atom::C_GLU_G)
    .value("C_GLU_D", atom::C_GLU_D)
    .value("O_GLU_E", atom::O_GLU_E)

    .value("O_ASN", atom::O_ASN)
    .value("N_ASN", atom::N_ASN)
    .value("C_ASN", atom::C_ASN)
    .value("C_ASN_A", atom::C_ASN_A)
    .value("C_ASN_B", atom::C_ASN_B)
    .value("C_ASN_G", atom::C_ASN_G)
    .value("O_ASN_D", atom::O_ASN_D)
    .value("N_ASN_D", atom::N_ASN_D)

    .value("O_GLN", atom::O_GLN)
    .value("N_GLN", atom::N_GLN)
    .value("C_GLN", atom::C_GLN)
    .value("C_GLN_A", atom::C_GLN_A)
    .value("C_GLN_B", atom::C_GLN_B)
    .value("C_GLN_G", atom::C_GLN_G)
    .value("C_GLN_D", atom::C_GLN_D)
    .value("O_GLN_E", atom::O_GLN_E)
    .value("N_GLN_E", atom::N_GLN_E)   

    .value("O_LYS", atom::O_LYS)
    .value("N_LYS", atom::N_LYS)
    .value("C_LYS", atom::C_LYS)
    .value("C_LYS_A", atom::C_LYS_A)
    .value("C_LYS_B", atom::C_LYS_B)
    .value("C_LYS_G", atom::C_LYS_G)
    .value("C_LYS_D", atom::C_LYS_D)
    .value("C_LYS_E", atom::C_LYS_E)
    .value("N_LYS_Z", atom::N_LYS_Z)     

    .value("O_ARG", atom::O_ARG)
    .value("N_ARG", atom::N_ARG)
    .value("C_ARG", atom::C_ARG)
    .value("C_ARG_A", atom::C_ARG_A)
    .value("C_ARG_B", atom::C_ARG_B)
    .value("C_ARG_G", atom::C_ARG_G)
    .value("C_ARG_D", atom::C_ARG_D)
    .value("N_ARG_E", atom::N_ARG_E)
    .value("C_ARG_Z", atom::C_ARG_Z)    
    .value("N_ARG_H", atom::N_ARG_H)

    .value("O_TYR", atom::O_TYR)
    .value("N_TYR", atom::N_TYR)
    .value("C_TYR", atom::C_TYR)
    .value("C_TYR_A", atom::C_TYR_A)
    .value("C_TYR_B", atom::C_TYR_B)
    .value("C_TYR_G", atom::C_TYR_G)
    .value("C_TYR_D", atom::C_TYR_D)
    .value("C_TYR_E", atom::C_TYR_E)
    .value("C_TYR_Z", atom::C_TYR_Z)
    .value("O_TYR_H", atom::O_TYR_H)

    .value("O_PHE", atom::O_PHE)
    .value("N_PHE", atom::N_PHE)
    .value("C_PHE", atom::C_PHE)
    .value("C_PHE_A", atom::C_PHE_A)
    .value("C_PHE_B", atom::C_PHE_B)
    .value("C_PHE_G", atom::C_PHE_G)
    .value("C_PHE_D", atom::C_PHE_D)
    .value("C_PHE_E", atom::C_PHE_E)
    .value("C_PHE_Z", atom::C_PHE_Z)

    .value("O_HIS", atom::O_HIS)
    .value("N_HIS", atom::N_HIS)
    .value("C_HIS", atom::C_HIS)
    .value("C_HIS_A", atom::C_HIS_A)
    .value("C_HIS_B", atom::C_HIS_B)
    .value("C_HIS_G", atom::C_HIS_G) 
    .value("N_HIS_D1", atom::N_HIS_D1)
    .value("C_HIS_D2", atom::C_HIS_D2)
    .value("C_HIS_E1", atom::C_HIS_E1)
    .value("N_HIS_E2", atom::N_HIS_E2)

    .value("O_TRP", atom::O_TRP)
    .value("N_TRP", atom::N_TRP)
    .value("C_TRP", atom::C_TRP)
    .value("C_TRP_A", atom::C_TRP_A)
    .value("C_TRP_B", atom::C_TRP_B)
    .value("C_TRP_G", atom::C_TRP_G) 
    .value("C_TRP_D1", atom::C_TRP_D1)
    .value("C_TRP_D2", atom::C_TRP_D2)
    .value("N_TRP_E1", atom::N_TRP_E1)
    .value("C_TRP_E2", atom::C_TRP_E2) 
    .value("C_TRP_E3", atom::C_TRP_E3) 
    .value("C_TRP_Z2", atom::C_TRP_Z2)
    .value("C_TRP_Z3", atom::C_TRP_Z3)     
    .value("C_TRP_H2", atom::C_TRP_H2)
    
    .value("O_PRO", atom::O_PRO)
    .value("N_PRO", atom::N_PRO)
    .value("C_PRO", atom::C_PRO)
    .value("C_PRO_A", atom::C_PRO_A)
    .value("C_PRO_B", atom::C_PRO_B)
    .value("C_PRO_G", atom::C_PRO_G)
    .value("C_PRO_D", atom::C_PRO_D)
  ;

  def("GetAtomTypeByName",&GetAtomTypeByName,(arg("amino_acid"),arg("atom_name")));

  class_<impl::InteractionOpts>("InteractionOpts", no_init)
    .def_readonly("lower_cutoff", &impl::InteractionOpts::lower_cutoff)
    .def_readonly("upper_cutoff", &impl::InteractionOpts::upper_cutoff)
    .def_readonly("number_of_bins", &impl::InteractionOpts::number_of_bins)
    .def_readonly("sequence_sep", &impl::InteractionOpts::sequence_sep)
    .def_readonly("sigma", &impl::InteractionOpts::sigma)
  ;

  class_<InteractionStatistic, bases<StatisticBase> > ("InteractionStatistic", init<Real, Real, int, int>())
    .def("Load", &InteractionStatistic::Load, (arg("filename"))).staticmethod("Load")
    .def("Save", &InteractionStatistic::Save, (arg("filename")))
    .def("Extract", &InteractionStatistic::Extract, (arg("target"),arg("env"),arg("weight")=1.0))
    .def("GetTotalCount", &InteractionStatistic::GetTotalCount)
    .def("GetCount", &GetCountatomdist, (arg("atom_a"), arg("atom_b"), arg("dist_bin")))
    .def("GetCount", &GetCountatom, (arg("atom_a"), arg("atom_b")))
    .def("GetCount", &GetCountdist, (arg("dist_bin")))
    .def("GetOptions", &InteractionStatistic::GetOpts)
  ;

  class_<InteractionPotential, bases<PotentialBase> >("InteractionPotential", no_init) //no_init=>forces to use create function
    .def("Load", &InteractionPotential::Load, (arg("filename"))).staticmethod("Load")
    .def("Save", &InteractionPotential::Save, (arg("filename")))
    .def("Create", &WrapHybridCreate, (arg("interaction_statistic_one"),arg("interaction_statistic_two"),arg("sigma")=0.02,arg("reference_state")="classic",arg("max_energy")=8.0))
    .def("Create", &WrapCreate, (arg("interaction_statistic"), arg("sigma")=0.02, arg("reference_state")="classic")).staticmethod("Create")    
    .def("SetEnvironment", &InteractionPotential::SetEnvironment, (arg("environment")))
    .def("GetEnergy", &GetEnergyatomdist, (arg("atom_a"),arg("atom_b"), arg("dist")))
    .def("GetEnergy", &GetEnergyresidueenv, (arg("target_residue"), arg("env"),arg("normalize")=true))
    .def("GetEnergy", &GetEnergyres, (arg("target_residue"),arg("normalize")=true))
    .def("GetEnergies", &WrapGetEnergies, (arg("target"), arg("environment"),arg("normalize")=true))
    .def("GetTotalEnergy", &InteractionPotential::GetTotalEnergy, (arg("target_view"), arg("environment_view"),arg("normalize")=true))
    .def("GetCounts", &InteractionPotential::GetCounts)
    .def("GetOptions", &InteractionPotential::GetOpts)
  ;



  register_ptr_to_python<InteractionStatisticPtr>(); 
  register_ptr_to_python<InteractionPotentialPtr>();
  
}

