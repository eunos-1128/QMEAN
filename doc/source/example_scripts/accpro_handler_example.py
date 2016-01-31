#Let's import the required module
from qmean import predicted_sequence_features
#we also need dssp, make sure an executable lies in your PATH!!
from ost.bindings import dssp


#let's define the accpro data for crambin
data = dict()
data["seq"] =  "TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN"
data["acc"] =  "bbbbbbeeeeeebbbbbebebbebebebbebbbbbbbbbeeeeebb"

#and generate an ACCPROHandler
handler = predicted_sequence_features.ACCPROHandler(data)

#let's load the crystal structure and assign the relative solvent
#acccessibility as estimated by dssp
prot = io.LoadPDB("1crn",remote=True)
dssp.AssignDSSP(prot, extract_burial_status=True)

#Let's use the handler for the full first chain, but also for a subset of it!
subset = prot.Select("rnum>5 and rnum<20")


#lets print out some stuff

print "Solvent accessibility for the full first chain:"
print handler.GetACCPROAccessibility(prot.chains[0])

print "Solvent accessibility for a subset:"
print handler.GetACCPROAccessibility(subset.chains[0])

print "And here are the according ACCAgreement scores"
print handler.GetACCAgreementFromChain(prot.chains[0], dssp_assigned=True)
print handler.GetACCAgreementFromChain(subset.chains[0], dssp_assigned=True)

