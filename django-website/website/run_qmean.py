import os, sys
from qmean import mqa_result
from qmean import mqa_result_membrane
from qmean.predicted_sequence_features import PSIPREDHandler
from qmean.predicted_sequence_features import ACCPROHandler
from ost.io import LoadPDB
from ost.io import SavePDB
from sm import config

USAGE = "usage: sm run_qmean.py <INPUT_DIR> <OUTPUT_DIR>"

if len(sys.argv) < 3:
  print USAGE
  sys.exit(-1)

in_dir = sys.argv[1]
out_dir = sys.argv[2]

#we only use the SGE module for submission, the run_qmean script takes
#now the responsibility over the status file
status_path = os.path.abspath(os.path.join(os.path.join(in_dir,os.pardir),"status"))
status_jobid = open(status_path,'r').readlines()[0].split()[0]

#let's change status to running
status_file = open(status_path,'w')
status_file.write(status_jobid)
status_file.write(" ")
status_file.write("RUNNING")
status_file.close()

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
                                dssp_path = config.DSSP_BIN)

  SavePDB(model,os.path.join(out_path,"model.pdb"))


#let's change the status to completed
status_file = open(status_path,'w')
status_file.write(status_jobid)
status_file.write(" ")
status_file.write("COMPLETED")
status_file.close()
