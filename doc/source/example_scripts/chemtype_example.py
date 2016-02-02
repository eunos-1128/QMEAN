#Load the required modules
from qmean import SSAgreement
from qmean import GetAtomTypeByName

#Let's check out all ChemTypes for a Phenylalanine...
atom_names = ["N","CA","CB","O","CG","CD1","CD2","CE1","CE2","CZ"]

for an in atom_names:
    chem_type = GetAtomTypeByName(ost.conop.PHE,an)
    print "ChemType for",an, ":\t", chem_type

#CD1 and CD2 have the same ChemType! The same is true for CE1 and CE2.
#The reason for that is that they cannot be distinguished, so QMEAN
#treats them equivalently.
