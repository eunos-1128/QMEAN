from ost import mol
from ost import seq
from ost import conop

def CheckClashes(ent):
  clashing_dist=mol.alg.DefaultClashingDistances()
  filtered_structure=mol.alg.FilterClashes(ent, clashing_dist, always_remove_bb=True)

  ret_list=list()

  for r in ent.residues:
    if filtered_structure.ViewForResidue(r.handle).IsValid():
      ret_list.append(0)
    else:
      ret_list.append(1)

  return ret_list

def CheckStereoChemistry(ent, std_bonds=12, std_angles=12):
  #be aware, that the default std parameters are rather generous
  bond_stereo_chemical_param=mol.alg.DefaultBondStereoChemicalParams()
  angle_stereo_chemical_param=mol.alg.DefaultAngleStereoChemicalParams()

  filtered_structure=mol.alg.CheckStereoChemistry(ent, bond_stereo_chemical_param, angle_stereo_chemical_param, std_bonds, std_angles, always_remove_bb=True)

  ret_list=list()

  for r in ent.residues:
    if filtered_structure.ViewForResidue(r.handle).IsValid():
      ret_list.append(0)
    else:
      ret_list.append(1) 

  return ret_list

def CheckResidueCompleteness(ent):

  builder=conop.GetDefaultLib()

  ret_list=list()

  for r in ent.residues:
    if builder.IsResidueComplete(r.handle):
      ret_list.append(0)
    else:
      ret_list.append(1)

  return ret_list

def SeqresCheck(ent, seqres):
  try:
    seq.alg.AlignToSEQRES(seqres)
    return True
  except:
    return False
