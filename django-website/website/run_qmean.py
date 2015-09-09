import os, sys
from qmean import mqa_result
from qmean import mqa_result_membrane
from qmean.predicted_sequence_features import PSIPREDHandler
from qmean.predicted_sequence_features import ACCPROHandler
from ost.io import LoadPDB

USAGE = "usage: ost run_qmean.py <INPUT_DIR> <OUTPUT_DIR>"

if len(sys.argv) < 3:
  print USAGE
  sys.exit(-1)

in_dir = sys.argv[1]
out_dir = sys.argv[2]

#at this point we would read the json file with all required input

files = os.listdir(in_dir)

model_files = list()
for f in files:
  if ".pdb" in f:
    model_files.append(f)


#sequence stuff has to be handled here
psipred_handler = None
accpro_handler = None


for f in model_files:

  model = LoadPDB(os.path.join(in_dir,f)).Select("peptide=true")
  out_path = os.path.join(out_dir,f.split('.')[0])
  mqa_result.AssessModelQuality(model, 
                                output_dir = out_path, 
                                psipred = psipred_handler,
                                accpro = accpro_handler,
                                dssp_path = "/import/bc2/apps/dssp/2003-May-30-goolf-1.4.10/bin/dssp")

