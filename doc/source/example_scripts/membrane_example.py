import os
from ost import io, mol, geom

#We use a GPCR example... 
#Directly after loading we strip of all hydrogens 
#and everything, that doesn't look like a peptide
gpcr_path = os.path.join("example_data","3EML.pdb")
gpcr = io.LoadPDB(gpcr_path).Select("peptide=true and ele!=H")

mem_param = mol.alg.FindMembrane(gpcr)

#The goal of this example is to create a structure only consisting
#of transmembrane residues. We take the CA position of every residue
#as reference position and check, whether it's inside the membrane planes
membrane_center = mem_param.pos * mem_param.membrane_axis
membrane_plane = geom.Plane(membrane_center, mem_param.membrane_axis)
half_width = mem_param.width * 0.5

for r in gpcr.residues:
  ca = r.FindAtom("CA")
  if not ca.IsValid():
    continue
  #Whether a residue is in the membrane or not can be defined
  #using generic properties of the ResidueHandle
  ca_pos = ca.GetPos()
  if abs(geom.Distance(membrane_plane, ca_pos)) < half_width:
    r.SetIntProp("in_membrane",1)
  else:
    r.SetIntProp("in_membrane",0)

#The generic properties can then be used to select the transmembrane part
transmembrane_part = gpcr.Select("grin_membrane:0=1")

#We can finally save down the transmembrane part!
io.SavePDB(gpcr,"gpcr.pdb")
io.SavePDB(transmembrane_part, "gpcr_transmembrane_part.pdb")
