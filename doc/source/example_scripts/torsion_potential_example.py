#Example code to demonstrate the usage of setting up a torsion
#potential with customized group identifiers and dihedral angle binning

import os
from ost import io
from qmean import TorsionStatistic, TorsionPotential


#Let's define some groups for three consecutive amino acids
group_identifier = list()

#If the central amino acid is Glycine, thats the group that gets chosen, 
#no matter what the neighbours are
group_identifier.append("all-GLY-all")

#Same is true for Proline
group_identifier.append("all-PRO-all")

#Also pre-Prolines get handled seperately but with lower priority
group_identifier.append("all-all-PRO")

#The last exception
group_identifier.append("all-VAL,ILE-all")

#Everything else gets threated by this group
group_identifier.append("all-all-all")

#Create a statistics object with the defined groups and
#bin sizes of 20 degrees for the central phi/psi
#torsion angles. One bin for the flanking amino acids
#basically neglects their torsion angles and only takes
#into account their identity.
bins = [1,1,18,18,1,1]
stat = TorsionStatistic(group_identifier, bins)

#Load some structure from the PDB and train the stats with it
crambin_path = os.path.join("example_data","1CRN.pdb")
crambin = io.LoadPDB(crambin_path).Select("cname=A and peptide=true")
stat.Extract(crambin)

#Create a super undersaturated potential
pot = TorsionPotential.Create(stat)

#Save the statistic and potential objects
stat.Save("torsion_stat.dat")
pot.Save("torsion_pot.dat")
