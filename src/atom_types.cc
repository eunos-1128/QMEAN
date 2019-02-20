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

#include <boost/version.hpp>
#if BOOST_VERSION>=103800
#define BOOST_SPIRIT_USE_OLD_NAMESPACE
#include <boost/spirit/include/classic_symbols.hpp>
#else
#include <boost/spirit/symbols.hpp>
#endif

#include <map>
#include "atom_types.hh"


using namespace ost::conop;

namespace qmean {

using namespace boost::spirit;


// let's have a little bit of syntactic sugar to define the atom name to
// atom type mapping.
class ResidueMap : public std::map<AminoAcid, atom::ChemType> {
public:
  ResidueMap(AminoAcid amino_acid, atom::ChemType type) {
    this->operator()(amino_acid, type);
  }
  ResidueMap& operator()(AminoAcid amino_acid, atom::ChemType type) {
    this->insert(std::make_pair(amino_acid, type));
    return *this;
  }
  atom::ChemType Get(AminoAcid amino_acid) {
    assert(amino_acid!=XXX);
    std::map<AminoAcid, atom::ChemType>::const_iterator i;
    i=this->find(amino_acid);
    if (i!=this->end())
      return i->second;
    return atom::UNKNOWN;
  }
};
typedef ResidueMap R;


struct AtomNames : public symbols<R> {
  AtomNames() {
    add
    ("C",   R(ALA, atom::C_ALA)
             (ARG, atom::C_ARG)
             (ASN, atom::C_ASN)
             (ASP, atom::C_ASP)
             (GLN, atom::C_GLN)
             (GLU, atom::C_GLU)
             (LYS, atom::C_LYS)
             (SER, atom::C_SER)
             (CYS, atom::C_CYS)
             (MET, atom::C_MET)
             (TRP, atom::C_TRP)
             (TYR, atom::C_TYR)
             (THR, atom::C_THR)
             (VAL, atom::C_VAL)
             (ILE, atom::C_ILE)
             (LEU, atom::C_LEU)
             (GLY, atom::C_GLY)
             (PRO, atom::C_PRO)
             (HIS, atom::C_HIS)
             (PHE, atom::C_PHE))
    ("N",   R(ALA, atom::N_ALA)
             (ARG, atom::N_ARG)
             (ASN, atom::N_ASN)
             (ASP, atom::N_ASP)
             (GLN, atom::N_GLN)
             (GLU, atom::N_GLU)
             (LYS, atom::N_LYS)
             (SER, atom::N_SER)
             (CYS, atom::N_CYS)
             (MET, atom::N_MET)
             (TRP, atom::N_TRP)
             (TYR, atom::N_TYR)
             (THR, atom::N_THR)
             (VAL, atom::N_VAL)
             (ILE, atom::N_ILE)
             (LEU, atom::N_LEU)
             (GLY, atom::N_GLY)
             (PRO, atom::N_PRO)
             (HIS, atom::N_HIS)
             (PHE, atom::N_PHE))
    ("O",   R(ALA, atom::O_ALA)
             (ARG, atom::O_ARG)
             (ASN, atom::O_ASN)
             (ASP, atom::O_ASP)
             (GLN, atom::O_GLN)
             (GLU, atom::O_GLU)
             (LYS, atom::O_LYS)
             (SER, atom::O_SER)
             (CYS, atom::O_CYS)
             (MET, atom::O_MET)
             (TRP, atom::O_TRP)
             (TYR, atom::O_TYR)
             (THR, atom::O_THR)
             (VAL, atom::O_VAL)
             (ILE, atom::O_ILE)
             (LEU, atom::O_LEU)
             (GLY, atom::O_GLY)
             (PRO, atom::O_PRO)
             (HIS, atom::O_HIS)
             (PHE, atom::O_PHE))
    ("CA",  R(ALA, atom::C_ALA_A)
             (ARG, atom::C_ARG_A)
             (ASN, atom::C_ASN_A)
             (ASP, atom::C_ASP_A)
             (GLN, atom::C_GLN_A)
             (GLU, atom::C_GLU_A)
             (LYS, atom::C_LYS_A)
             (SER, atom::C_SER_A)
             (CYS, atom::C_CYS_A)
             (MET, atom::C_MET_A)
             (TRP, atom::C_TRP_A)
             (TYR, atom::C_TYR_A)
             (THR, atom::C_THR_A)
             (VAL, atom::C_VAL_A)
             (GLY, atom::C_GLY_A)
             (ILE, atom::C_ILE_A)
             (LEU, atom::C_LEU_A)
             (PRO, atom::C_PRO_A)
             (HIS, atom::C_HIS_A)
             (PHE, atom::C_PHE_A))
    ("CB",  R(ALA, atom::C_ALA_B)
             (ARG, atom::C_ARG_B)
             (ASN, atom::C_ASN_B)
             (ASP, atom::C_ASP_B)
             (GLN, atom::C_GLN_B)
             (GLU, atom::C_GLU_B)
             (LYS, atom::C_LYS_B)
             (SER, atom::C_SER_B)
             (CYS, atom::C_CYS_B)
             (MET, atom::C_MET_B)
             (TRP, atom::C_TRP_B)
             (TYR, atom::C_TYR_B)
             (THR, atom::C_THR_B)
             (VAL, atom::C_VAL_B)
             (ILE, atom::C_ILE_B)
             (LEU, atom::C_LEU_B)
             (PRO, atom::C_PRO_B)
             (HIS, atom::C_HIS_B)
             (PHE, atom::C_PHE_B))
    ("CG",  R(ARG, atom::C_ARG_G)
             (ASN, atom::C_ASN_G)
             (ASP, atom::C_ASP_G)
             (GLN, atom::C_GLN_G)
             (GLU, atom::C_GLU_G)
             (LYS, atom::C_LYS_G)
             (MET, atom::C_MET_G)
             (TRP, atom::C_TRP_G)
             (LEU, atom::C_LEU_G)
             (PRO, atom::C_PRO_G)
             (PHE, atom::C_PHE_G)
             (TYR, atom::C_TYR_G)
             (HIS, atom::C_HIS_G))
    ("CG1", R(ILE, atom::C_ILE_G1)
             (VAL, atom::C_VAL_G))
    ("CG2", R(ILE, atom::C_ILE_G2)
             (THR, atom::C_THR_G2)
             (VAL, atom::C_VAL_G))
    ("CD",  R(GLU, atom::C_GLU_D)
             (ILE, atom::C_GLN_D)
             (LYS, atom::C_LYS_D)
             (PRO, atom::C_PRO_D)
             (GLN, atom::C_GLN_D)
             (GLU, atom::C_GLU_D)
             (ARG, atom::C_ARG_D))
    ("CD1", R(TRP, atom::C_TRP_D1)
             (TYR, atom::C_TYR_D)
             (PHE, atom::C_PHE_D)
             (LEU, atom::C_LEU_D)
             (ILE, atom::C_ILE_D1))
    ("CD2", R(TRP, atom::C_TRP_D2)
             (TYR, atom::C_TYR_D)
             (PHE, atom::C_PHE_D)
             (LEU, atom::C_LEU_D)
             (HIS, atom::C_HIS_D2))
    ("CE",  R(LYS, atom::C_LYS_E)
             (MET, atom::C_MET_E))
    ("CE1", R(PHE, atom::C_PHE_E)
             (TYR, atom::C_TYR_E)
             (HIS, atom::C_HIS_E1))
    ("CE2", R(PHE, atom::C_PHE_E)
             (TYR, atom::C_TYR_E)
             (TRP, atom::C_TRP_E2))
    ("CE3", R(TRP, atom::C_TRP_E3))
    ("CH2", R(TRP, atom::C_TRP_H2))
    ("CZ",  R(ARG, atom::C_ARG_Z)
             (PHE, atom::C_PHE_Z)
             (TYR, atom::C_TYR_Z))
    ("CZ2", R(TRP, atom::C_TRP_Z2))
    ("CZ3", R(TRP, atom::C_TRP_Z3))
    ("OG",  R(SER, atom::O_SER_G))
    ("OG1", R(THR, atom::O_THR_G1))
    ("SG",  R(CYS, atom::S_CYS_G))
    ("SD",  R(MET, atom::S_MET_D))
    ("OD1", R(ASP, atom::O_ASP_D)
             (ASN, atom::O_ASN_D))
    ("OD2", R(ASP, atom::O_ASP_D))
    ("OE1", R(GLU, atom::O_GLU_E)
             (GLN, atom::O_GLN_E))
    ("OE2", R(GLU, atom::O_GLU_E))
    ("ND1", R(HIS, atom::N_HIS_D1))
    ("ND2", R(ASN, atom::N_ASN_D))
    ("NE",  R(ARG, atom::N_ARG_E))
    ("NE1", R(ARG, atom::N_ARG_E)
             (TRP, atom::N_TRP_E1))
    ("NE2", R(GLN, atom::N_GLN_E)
             (HIS, atom::N_HIS_E2))
    ("NH1", R(ARG, atom::N_ARG_H))
    ("NH2", R(ARG, atom::N_ARG_H))
    ("NZ",  R(LYS, atom::N_LYS_Z))
    ("OH",  R(TYR, atom::O_TYR_H))
     ;
  }
};

atom::ChemType GetAtomTypeByName(AminoAcid amino_acid, const String& aname) {

  static AtomNames aa_name_symbols;
  R* m=find(aa_name_symbols, aname.c_str());
  if (m) {
    atom::ChemType t=m->Get(amino_acid);
    return t;
  }
  return atom::UNKNOWN;

}
}
