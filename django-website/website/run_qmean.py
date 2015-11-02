import os, sys, json, shutil, fcntl
try:
  import hashlib
except ImportError:
  import md5 as hashlib

import matplotlib
matplotlib.use("Agg")
from qmean import mqa_result
from qmean import mqa_result_membrane
from qmean.predicted_sequence_features import PSIPREDHandler
from qmean.predicted_sequence_features import ACCPROHandler
from ost.io import LoadPDB
from ost.io import SavePDB
from ost.io import LoadSequence
from ost.io import SaveSequence
from ost.bindings import hhblits
from ost.seq import CreateSequence
from sm import config
from sm.smtl import structanno


def CopyFromCache(cache_dir,hash,path):

  cache_file = os.path.join(cache_dir,hash)
  lock = os.path.join(cache_dir,hash+".lock")

  if not os.path.exists(lock):
    if os.path.exists(cache_file):
      try:
        shutil.copy(cache_file, path)
        return True
      except:
        pass
  return False


def CopyToCache(cache_dir,hash,path):

  cache_file = os.path.join(cache_dir,hash)
  lock = os.path.join(cache_dir,hash+".lock")

  if not os.path.exists(cache_file) and not os.path.exists(lock):
    lock_fd = open(lock, 'w')
    try:
      fcntl.lockf(lock_fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except IOError:
      pass
    else:
      try:
        shutil.copy(path, cache_file)
      except:
        pass
      fcntl.lockf(lock_fd, fcntl.LOCK_UN)
      lock_fd.close()
      try:
        os.remove(lock)
      except:
        pass


def GetSequencePredictions(sequences, hhblits_cache, accpro_cache, out_dir):

  a3m_files = list()
  acc_files = list()

  sequence_hashes = list()
  for s in sequences:
    sequence_hashes.append(hashlib.md5(s.GetString()).hexdigest())  

  #fire HHBlits
  for s,h in zip(sequences,sequence_hashes):
    a3m_file = os.path.join(out_dir,s.GetName()+".a3m")
    if CopyFromCache(hhblits_cache, h, a3m_file):
      a3m_files.append(a3m_file)
      continue 
    hh = hhblits.HHblits(s, hhsuite_root = config.HHSUITE_ROOT_DIR,
                         working_dir=out_dir)

    query_msa = hh.BuildQueryMSA(nrdb = config.HHBLITS_NRDB_PREFIX, 
                                 cpu = config.HHBLITS_CPU)
    shutil.copy(query_msa, a3m_file)
    a3m_files.append(a3m_file)
    CopyToCache(hhblits_cache,h,a3m_file)

  #read in the profiles from the a3m files
  profiles = list()
  for f in a3m_files:
    f_h = open(f,'r')
    content = f_h.readlines()
    f_h.close()
    profiles.append(hhblits.ParseA3M(content))

  #fire ACCPRO...
  for s,h,p in zip(sequences,sequence_hashes,profiles):
    accpro_file = os.path.join(out_dir,s.GetName()+".acc")
    if CopyFromCache(accpro_cache,h,accpro_file):
      acc_files.append(accpro_file)
      continue
    target_fasta_path = os.path.join(out_dir,"accpro_target.fasta") 
    SaveSequence(s, target_fasta_path)
    aln_filename = structanno._PrepareSS_ACCPro(p['msa'])
    _, accpro = structanno._GetSS_ACCproAnnotation(os.path.join(out_dir,"accpro_target.fasta"),
                                                 aln_filename,
                                                 s.GetLength())
    os.remove(target_fasta_path)
    accpro = "".join(accpro)
    seq_to_save = CreateSequence("accpro",accpro)
    SaveSequence(seq_to_save,accpro_file,format="fasta")
    acc_files.append(accpro_file)
    CopyToCache(accpro_cache,h,accpro_file)

  #lets gather the final predictions
  psipred_pred = list()
  psipred_conf = list()
  accpro_pred = list()

  for i,p in enumerate(profiles):
    psipred_pred.append("".join([item for item in p["ss_pred"]]))
    psipred_conf.append("".join([str(item) for item in p["ss_conf"]]))
    acc_seq = LoadSequence(acc_files[i], format = "fasta")
    accpro_pred.append(acc_seq.GetString())

  return (psipred_pred, psipred_conf, accpro_pred, a3m_files, acc_files)


def GetDistanceConstraints(sequences, a3m_files):
  return list()

def MatchSeqPredictionsWithModel(sequences, model, psipred_pred,
                                 psipred_conf, accpro_pred):
  

  psipred_handler = list()
  accpro_handler = list()

  #The input data has been checked... so we assume that the sequences
  #match with the chains given one of the three following scenarios:

  #If there is no sequence, there won't be any sequence prediction...
  if len(sequences) == 0:
    return None, None
  #if there is one sequence, the model either has one chain or
  #the model is a homo-oligomer
  elif len(sequences) == 1:
    data = dict()
    data["seq"] = sequences[0].GetString()
    data["ss"] = psipred_pred[0]
    data["conf"] = psipred_conf[0]
    data["acc"] = accpro_pred[0]
    for c in model.chains:
      pp = PSIPREDHandler(data)
      acc = ACCPROHandler(data)
      psipred_handler.append(pp)
      accpro_handler.append(acc)
  #if there is more than one sequence, we match the chain/sequence names
  else:
    for c in model.chains:
      c_name = c.GetName()
      found_seq = False
      for i,s in enumerate(sequences):
        if s.GetName() == c_name:
          data = dict()
          data["seq"] = s.GetString()
          data["ss"] = psipred_pred[i]
          data["conf"] = psipred_conf[i]
          data["acc"] = accpro_pred[i]
          pp = PSIPREDHandler(data)
          acc = ACCPROHandler(data)
          psipred_handler.append(pp)
          accpro_handler.append(acc)
          found_seq = True
          break
      if not found_seq:
        raise RuntimeError("Could not find matching sequence for chain of name "+c_name)

  return (psipred_handler, accpro_handler)


def MatchDCWithModel(sequences,model,dc_list):
  return None


USAGE = "usage: sm run_qmean.py <INPUT_DIR> <OUTPUT_DIR> <CACHE_DIR>"

if len(sys.argv) < 4:
  print USAGE
  sys.exit(-1)

in_dir = sys.argv[1]
out_dir = sys.argv[2]
cache_dir = sys.argv[3]

if not os.path.exists(in_dir):
  raise ValueError("Provided input directory does not exist!")
if not os.path.exists(out_dir):
  raise ValueError("Provided output directory does not exist!")
if not os.path.exists(cache_dir):
  raise ValueError("Provided cache directory does not exist!")

#We build caches for the hhblits profiles and the accpro results...
hhblits_cache = os.path.join(cache_dir,"hhblits-cache")
accpro_cache = os.path.join(cache_dir,"accpro-cache")

if not os.path.exists(hhblits_cache):
  os.makedirs(hhblits_cache)
if not os.path.exists(accpro_cache):
  os.makedirs(accpro_cache)

#we only use the SGE module for submission, the run_qmean script takes
#now the responsibility over the status file. Since the whole thing is sequential
#and only one job per project is running we don't care for concurrent write issues
status_path = os.path.abspath(os.path.join(os.path.join(in_dir,os.pardir),"status"))
status_jobid = open(status_path,'r').readlines()[0].split()[0]

#let's change status to running
status_file = open(status_path,'w')
msg = status_jobid + " " + "RUNNING"
status_file.write(msg)
status_file.close()

#lets read the json stuff to get all required information about the project
infile = open(os.path.join(in_dir,"project.json"))
project_data = json.load(infile)
infile.close()

use_qmeandisco = project_data["options"]["qmeandisco"]
use_qmeanbrane = project_data["options"]["qmeanbrane"]

#lets extract the names of the models
model_files = list()
for item in project_data["models"]:
  model_files.append(str(item["modelid"])+".pdb")

#let's extract all sequences
sequences = list()
for item in project_data["sequences"]:
  sequences.append(ost.seq.CreateSequence(str(item["name"]),str(item["sequence"])))

#Let's get all stuff we predict only with the sequence
psipred_pred, psipred_conf, accpro_pred, a3m_files, accpro_files = GetSequencePredictions(sequences, 
                                                                                          hhblits_cache, 
                                                                                          accpro_cache, 
                                                                                          out_dir)

db_list = list()
if use_qmeandisco:
  dc_list = GetDistanceConstraints(sequences, a3m_files)


for f in model_files:

  model = LoadPDB(os.path.join(in_dir,f)).Select("peptide=true")

  psipred_handler, accpro_handler = MatchSeqPredictionsWithModel(sequences, 
                                                                 model, 
                                                                 psipred_pred,
                                                                 psipred_conf, 
                                                                 accpro_pred)

  if use_qmeandisco:
    dc = MatchDCWithModel(sequences, model, dc_list)
  else:
    dc = None

  out_path = os.path.join(out_dir,f.split('.')[0])
  mqa_result.AssessModelQuality(model, 
                                output_dir = out_path, 
                                psipred = psipred_handler,
                                accpro = accpro_handler,
                                dssp_path = config.DSSP_BIN)

  SavePDB(model,os.path.join(out_path,"model.pdb"))


#let's change the status to completed
status_file = open(status_path,'w')
msg = status_jobid + " " + "COMPLETED"
status_file.write(msg)
status_file.close()

