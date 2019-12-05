import os
from ost import io, mol
from qmean import PSIPREDHandler

#let's define the psipred data for crambin
data = dict()
data["seq"] =  "TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN"
data["ss"] =   "CCCCCCHHHHHHHHHEECCCCCHHHHHHHCCCEECCCCCCCCCCCC"
data["conf"] = "9657578887665410368988878977449676559978188898"

#and generate a PSIPREDHandler
handler = PSIPREDHandler(data)

#let's load the crystal structure and assign the secondary structure
#as estimated by the dssp algorithm
prot_path = os.path.join("example_data","1CRN.pdb")
prot = io.LoadPDB(prot_path)
mol.alg.AssignSecStruct(prot)

#Let's use the handler for the full first chain, but also for a subset of it!
subset = prot.Select("rnum>5 and rnum<20")


#lets print out some stuff

print "Psipred prediction and confidence for the full first chain:"
print handler.GetPSIPREDSS(prot.chains[0])
print handler.GetPSIPREDConf(prot.chains[0])

print "Psipred prediction and confidence for a subset:"
print handler.GetPSIPREDSS(subset.chains[0])
print handler.GetPSIPREDConf(subset.chains[0])

print "And here are the according SSAgreement scores"
print handler.GetSSAgreementFromChain(prot.chains[0])
print handler.GetSSAgreementFromChain(subset.chains[0])

