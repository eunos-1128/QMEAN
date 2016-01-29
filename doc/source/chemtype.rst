The Atom ChemType
================================================================================


Some functionality in QMEAN requires to have a definition of all chemically
distinguishable heavy atoms in protein structures. This is achieved by defining
the ChemType enum. You can either directly access them or generate them 
using following function:

.. method:: GetAtomTypeByName(aa, aname)

  :param aa:            Amino Acid the atom belongs to
  :param aname:         Default PDB name of the atom

  :type aa:             :class:`ost.conop.AminoAcid`
  :type aname:          :class:`str`

  :returns:             The enum type defining the according heavy atom,
                        ChemType.UNKNOWN if not existing.
  :rtype:               :ref:`ChemType <ChemType>`         



.. _ChemType:

* The ChemType Enum:

  * ChemType.O_GLY
  * ChemType.N_GLY
  * ChemType.C_GLY    
  * ChemType.C_GLY_A
  * ChemType.O_ALA 
  * ChemType.N_ALA 
  * ChemType.C_ALA     
  * ChemType.C_ALA_A
  * ChemType.C_ALA_B
  * ChemType.O_VAL
  * ChemType.N_VAL
  * ChemType.C_VAL    
  * ChemType.C_VAL_A
  * ChemType.C_VAL_B
  * ChemType.C_VAL_G
  * ChemType.O_LEU
  * ChemType.N_LEU
  * ChemType.C_LEU  
  * ChemType.C_LEU_A
  * ChemType.C_LEU_B
  * ChemType.C_LEU_G
  * ChemType.C_LEU_D  
  * ChemType.O_ILE
  * ChemType.N_ILE
  * ChemType.C_ILE
  * ChemType.C_ILE_A
  * ChemType.C_ILE_B
  * ChemType.C_ILE_G1
  * ChemType.C_ILE_G2
  * ChemType.C_ILE_D1
  * ChemType.O_THR
  * ChemType.N_THR
  * ChemType.C_THR
  * ChemType.C_THR_A
  * ChemType.C_THR_B
  * ChemType.O_THR_G1
  * ChemType.C_THR_G2
  * ChemType.O_SER
  * ChemType.N_SER
  * ChemType.C_SER
  * ChemType.C_SER_A
  * ChemType.C_SER_B
  * ChemType.O_SER_G
  * ChemType.O_CYS
  * ChemType.N_CYS
  * ChemType.C_CYS
  * ChemType.C_CYS_A
  * ChemType.C_CYS_B
  * ChemType.S_CYS_G    
  * ChemType.O_MET
  * ChemType.N_MET
  * ChemType.C_MET
  * ChemType.C_MET_A
  * ChemType.C_MET_B
  * ChemType.C_MET_G
  * ChemType.S_MET_D
  * ChemType.C_MET_E
  * ChemType.O_ASP
  * ChemType.N_ASP
  * ChemType.C_ASP
  * ChemType.C_ASP_A
  * ChemType.C_ASP_B
  * ChemType.C_ASP_G
  * ChemType.O_ASP_D
  * ChemType.O_GLU
  * ChemType.N_GLU
  * ChemType.C_GLU
  * ChemType.C_GLU_A
  * ChemType.C_GLU_B
  * ChemType.C_GLU_G
  * ChemType.C_GLU_D
  * ChemType.O_GLU_E
  * ChemType.O_ASN
  * ChemType.N_ASN
  * ChemType.C_ASN
  * ChemType.C_ASN_A
  * ChemType.C_ASN_B
  * ChemType.C_ASN_G
  * ChemType.O_ASN_D
  * ChemType.N_ASN_D
  * ChemType.O_GLN
  * ChemType.N_GLN
  * ChemType.C_GLN
  * ChemType.C_GLN_A
  * ChemType.C_GLN_B
  * ChemType.C_GLN_G
  * ChemType.C_GLN_D
  * ChemType.O_GLN_E
  * ChemType.N_GLN_E   
  * ChemType.O_LYS
  * ChemType.N_LYS
  * ChemType.C_LYS
  * ChemType.C_LYS_A
  * ChemType.C_LYS_B
  * ChemType.C_LYS_G
  * ChemType.C_LYS_D
  * ChemType.C_LYS_E
  * ChemType.N_LYS_Z     
  * ChemType.O_ARG
  * ChemType.N_ARG
  * ChemType.C_ARG
  * ChemType.C_ARG_A
  * ChemType.C_ARG_B
  * ChemType.C_ARG_G
  * ChemType.C_ARG_D
  * ChemType.N_ARG_E
  * ChemType.C_ARG_Z    
  * ChemType.N_ARG_H
  * ChemType.O_TYR
  * ChemType.N_TYR
  * ChemType.C_TYR
  * ChemType.C_TYR_A
  * ChemType.C_TYR_B
  * ChemType.C_TYR_G
  * ChemType.C_TYR_D
  * ChemType.C_TYR_E
  * ChemType.C_TYR_Z
  * ChemType.O_TYR_H
  * ChemType.O_PHE
  * ChemType.N_PHE
  * ChemType.C_PHE
  * ChemType.C_PHE_A
  * ChemType.C_PHE_B
  * ChemType.C_PHE_G
  * ChemType.C_PHE_D
  * ChemType.C_PHE_E
  * ChemType.C_PHE_Z
  * ChemType.O_HIS
  * ChemType.N_HIS
  * ChemType.C_HIS
  * ChemType.C_HIS_A
  * ChemType.C_HIS_B
  * ChemType.C_HIS_G 
  * ChemType.N_HIS_D1
  * ChemType.C_HIS_D2
  * ChemType.C_HIS_E1
  * ChemType.N_HIS_E2
  * ChemType.O_TRP
  * ChemType.N_TRP
  * ChemType.C_TRP
  * ChemType.C_TRP_A
  * ChemType.C_TRP_B
  * ChemType.C_TRP_G 
  * ChemType.C_TRP_D1
  * ChemType.C_TRP_D2
  * ChemType.N_TRP_E1
  * ChemType.C_TRP_E2 
  * ChemType.C_TRP_E3 
  * ChemType.C_TRP_Z2
  * ChemType.C_TRP_Z3     
  * ChemType.C_TRP_H2
  * ChemType.O_PRO
  * ChemType.N_PRO
  * ChemType.C_PRO
  * ChemType.C_PRO_A
  * ChemType.C_PRO_B
  * ChemType.C_PRO_G
  * ChemType.C_PRO_D
  * ChemType.UNKNOWN
















