# QMEAN:
# ----------
# 
# Copyright (C) 2018 University of Basel and the Authors.
# Authors: Gabriel Studer, Pascal Benkert, Marco Biasini
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

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
