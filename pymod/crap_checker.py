# Copyright (c) 2013-2018, SIB - Swiss Institute of Bioinformatics and
# Biozentrum - University of Basel
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
# http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

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
