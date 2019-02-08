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

#ifndef QMEAN_ATOM_TYPES_HH
#define QMEAN_ATOM_TYPES_HH

#include <qmean/module_config.hh>
#include <ost/conop/amino_acids.hh>

namespace qmean {

namespace atom {

  // list of chemically distinguishable atoms. 
  typedef enum {
    O_GLY = 0,
    N_GLY,
    C_GLY,    
    C_GLY_A,

    O_ALA,
    N_ALA,
    C_ALA,    
    C_ALA_A,
    C_ALA_B,

    O_VAL,
    N_VAL,
    C_VAL,    
    C_VAL_A,
    C_VAL_B,
    C_VAL_G,

    O_LEU,
    N_LEU,
    C_LEU,  
    C_LEU_A,
    C_LEU_B,
    C_LEU_G,
    C_LEU_D,  

    O_ILE,
    N_ILE,
    C_ILE,
    C_ILE_A,
    C_ILE_B,
    C_ILE_G1,
    C_ILE_G2,
    C_ILE_D1,

    O_THR,
    N_THR,
    C_THR,
    C_THR_A,
    C_THR_B,
    O_THR_G1,
    C_THR_G2,

    O_SER,
    N_SER,
    C_SER,
    C_SER_A,
    C_SER_B,
    O_SER_G,

    O_CYS,
    N_CYS,
    C_CYS,
    C_CYS_A,
    C_CYS_B,
    S_CYS_G,    

    O_MET,
    N_MET,
    C_MET,
    C_MET_A,
    C_MET_B,
    C_MET_G,
    S_MET_D,
    C_MET_E,

    O_ASP,
    N_ASP,
    C_ASP,
    C_ASP_A,
    C_ASP_B,
    C_ASP_G,
    O_ASP_D,

    O_GLU,
    N_GLU,
    C_GLU,
    C_GLU_A,
    C_GLU_B,
    C_GLU_G,
    C_GLU_D,
    O_GLU_E,

    O_ASN,
    N_ASN,
    C_ASN,
    C_ASN_A,
    C_ASN_B,
    C_ASN_G,
    O_ASN_D,
    N_ASN_D,

    O_GLN,
    N_GLN,
    C_GLN,
    C_GLN_A,
    C_GLN_B,
    C_GLN_G,
    C_GLN_D,
    O_GLN_E,
    N_GLN_E,   

    O_LYS,
    N_LYS,
    C_LYS,
    C_LYS_A,
    C_LYS_B,
    C_LYS_G,
    C_LYS_D,
    C_LYS_E,
    N_LYS_Z,     

    O_ARG,
    N_ARG,
    C_ARG,
    C_ARG_A,
    C_ARG_B,
    C_ARG_G,
    C_ARG_D,
    N_ARG_E,
    C_ARG_Z,    
    N_ARG_H,

    O_TYR,
    N_TYR,
    C_TYR,
    C_TYR_A,
    C_TYR_B,
    C_TYR_G,
    C_TYR_D, //D1 AND D2
    C_TYR_E, //E1 and E2
    C_TYR_Z,
    O_TYR_H,

    O_PHE,
    N_PHE,
    C_PHE,
    C_PHE_A,
    C_PHE_B,
    C_PHE_G,
    C_PHE_D, //D1 AND D2
    C_PHE_E, //E1 AND E2
    C_PHE_Z,

    O_HIS,
    N_HIS,
    C_HIS,
    C_HIS_A,
    C_HIS_B,
    C_HIS_G, 
    N_HIS_D1,
    C_HIS_D2,
    C_HIS_E1,
    N_HIS_E2,

    O_TRP,
    N_TRP,
    C_TRP,
    C_TRP_A,
    C_TRP_B,
    C_TRP_G, 
    C_TRP_D1,
    C_TRP_D2,
    N_TRP_E1,
    C_TRP_E2, 
    C_TRP_E3, 
    C_TRP_Z2,
    C_TRP_Z3,     
    C_TRP_H2,
    
    O_PRO,
    N_PRO,
    C_PRO,
    C_PRO_A,
    C_PRO_B,
    C_PRO_G,
    C_PRO_D,
    UNKNOWN
  } ChemType;
}
 
 
DLLEXPORT_QMEAN atom::ChemType GetAtomTypeByName(ost::conop::AminoAcid aa, 
                                                 const String& aname);

}

#endif

