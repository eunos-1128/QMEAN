import os
from ost.io import LoadPDB
from qmean import AssessMembraneModelQuality

#we load and assess some alternative models for the same target as they
#are used as a modelling example in the QMEANBrane publication
model_names = ["original_hhblits_alignment", 
               "shift_in_front_helix_four",
               "shift_into_middle",
               "shift_towards_cter"]

for mn in model_names:
  model_path = os.path.join("example_data",mn+".pdb")
  model = LoadPDB(model_path)
  AssessMembraneModelQuality(model,output_dir=mn)
