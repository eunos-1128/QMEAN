#Load QMEAN modules
from qmean import *
#import python random module
import random
import os

#Load crambin and extract CA positions
crambin_path = os.path.join("example_data","1CRN.pdb")
crambin = io.LoadPDB(crambin_path)
positions = [r.FindAtom("CA").GetPos() for r in crambin.residues]

#Generate smoother object with extracted positions
#and a std of 5.0
smoother = SphericalSmoother(positions, 5.0)

#Generate a random value for every residue in the structure
values = [random.random() for i in range(len(crambin.residues))]

#Smooth the random values
smoothed_values = smoother.Smooth(values)

print "initial values: ",values
print
print "smoothed values: ", smoothed_values
