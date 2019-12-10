import os
from ost import io, mol
from qmean import ACCPROHandler

#let's define the accpro data for crambin
data = dict()
data["seq"] =  "TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN"
data["acc"] =  "bbbbbbeeeeeebbbbbebebbebebebbebbbbbbbbbeeeeebb"

#and generate an ACCPROHandler
handler = ACCPROHandler(data)

#let's load the crystal structure and assign the relative solvent
#acccessibility
crambin_path = os.path.join("example_data","1CRN.pdb")
crambin = io.LoadPDB(crambin_path)
mol.alg.Accessibility(crambin, algorithm = mol.alg.DSSP)

#Let's use the handler for the full first chain, but also for a subset of it!
subset = crambin.Select("rnum>5 and rnum<20")

#lets print out some stuff
print("Solvent accessibility for the full first chain:")
print(handler.GetACCPROAccessibility(crambin.chains[0]))

print("Solvent accessibility for a subset:")
print(handler.GetACCPROAccessibility(subset.chains[0]))

print("And here are the according ACCAgreement scores")
print(handler.GetACCAgreementFromChain(crambin.chains[0]))
print(handler.GetACCAgreementFromChain(subset.chains[0]))
