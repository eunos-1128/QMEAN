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
from qmean import FindMembrane
from ost.io import LoadPDB
from ost.io import SavePDB
from ost.io import LoadSequence
from ost.io import SaveSequence
from ost.bindings import hhblits
from ost.bindings import naccess
from ost.bindings import msms
from ost.seq import CreateSequence
from ost.seq import CreateSequenceList
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

def CreateMembraneRepresentation(ent, mem_param = None, membrane_margin = 15, delta = 2.0):

  if mem_param == None:
    #let's calculate the transmembrane part of the structure
    membrane_finder_input = ent.Select("peptide=true and ele!=H")
    surf = msms.CalculateSurface(membrane_finder_input,radius=1.4)[0]
    naccess.CalculateSurfaceArea(membrane_finder_input)
    asa = list()
    for a in membrane_finder_input.atoms:
      if a.HasProp('asaAtom'):
        asa.append(a.GetFloatProp('asaAtom'))
      else:
        asa.append(0.0)
    mem_param = FindMembrane(mol.CreateEntityFromView(membrane_finder_input,False), surf, asa)

  membrane_axis = geom.Normalize(mem_param.membrane_axis)
  membrane_pos = mem_param.pos
  membrane_width = mem_param.width

  #define planes, that define the membrane
  top = membrane_pos * membrane_axis + membrane_axis * membrane_width/2
  bottom = membrane_pos * membrane_axis - membrane_axis * membrane_width/2
  plane_one = geom.Plane(top, membrane_axis)
  plane_two = geom.Plane(bottom, membrane_axis)

  #find residues close to those planes to find the actual membrane location
  for r in ent.residues:
    ca = r.FindAtom("CA")
    if not ca.IsValid():
      continue
    dist_one = abs(geom.Distance(plane_one,ca.GetPos()))
    dist_two = abs(geom.Distance(plane_two,ca.GetPos()))

    if dist_one < 3:
      r.SetIntProp("close_plane",1)
      continue

    if dist_two < 3:
      r.SetIntProp("close_plane",1)

  close_to_plane_selection = ent.Select("grclose_plane:0=1")
  membrane_center = close_to_plane_selection.GetGeometricCenter()

  #find the residue being most away from the center of the membrane
  max_dist_to_membrane_center = 0.0
  center_line = geom.Line3(membrane_center,membrane_center+membrane_axis)
  for r in close_to_plane_selection.residues:
    ca = r.FindAtom("CA")
    if not ca.IsValid():
      continue
    actual_dist = geom.Distance(center_line,ca.GetPos())
    max_dist_to_membrane_center = max(max_dist_to_membrane_center, actual_dist)


  #reassign the top and bottom positions, that have been only arbitrary points on the membrane planes
  top = geom.IntersectionPoint(center_line,plane_one)
  bottom = geom.IntersectionPoint(center_line,plane_two)


  #find a pair of perpendicular vectors, that are on the plane
  arbitrary_vec = geom.Vec3()
  if abs(membrane_axis[0]) <= abs(membrane_axis[1]) and abs(membrane_axis[0]) <= abs(membrane_axis[2]):
    arbitrary_vec = geom.Vec3(1,0,0)
  elif abs(membrane_axis[1]) <= abs(membrane_axis[0]) and abs(membrane_axis[1]) <= abs(membrane_axis[2]): 
    arbitrary_vec = geom.Vec3(0,1,0)
  else:
    arbitrary_vec = geom.Vec3(0,0,1)
  plane_x = geom.Normalize(geom.Cross(membrane_axis, arbitrary_vec))
  plane_y = geom.Normalize(geom.Cross(membrane_axis, plane_x))


  #Find the dummy atom positions
  radius = max_dist_to_membrane_center + membrane_margin
  num_sampling_points = int(round(2*radius/delta))
  mem_positions = list()

  #do plane one
  origin = top - delta * num_sampling_points / 2 * plane_x - delta * num_sampling_points / 2 * plane_y
  for i in range(num_sampling_points):
    for j in range(num_sampling_points):
      pos = origin + i * delta * plane_x + j * delta * plane_y
      if geom.Distance(pos,top) > radius:
        continue
      if len(ent.FindWithin(pos,4.0)) > 0:
        continue
      mem_positions.append(pos)

  #do plane two
  origin = bottom - delta * num_sampling_points / 2 * plane_x - delta * num_sampling_points / 2 * plane_y
  for i in range(num_sampling_points):
    for j in range(num_sampling_points):
      pos = origin + i * delta * plane_x + j * delta * plane_y
      if geom.Distance(pos,bottom) > radius:
        continue
      if len(ent.FindWithin(pos,4.0)) > 0:
        continue
      mem_positions.append(pos)

  if len(mem_positions) >= 9999:
    raise RuntimeError("Need more fancy implementation to get a dummy membrane entity, that can be saved down to disk")

  #fill positions into the entity
  mem_entity = mol.CreateEntity()
  ed = mem_entity.EditXCS()
  chain = ed.InsertChain("M")

  for p in mem_positions:
    res = ed.AppendResidue(chain,"MEM")
    ed.InsertAtom(res,"MEM",p)
      
  return mem_entity



def GetSequenceFeatures(sequences, hhblits_cache, accpro_cache, out_dir):

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
  accpro_failed = False
  for s,h,p in zip(sequences,sequence_hashes,profiles):
    try:
      accpro_file = os.path.join(out_dir,s.GetName()+".acc")
      if CopyFromCache(accpro_cache,h,accpro_file):
        acc_files.append(accpro_file)
        continue
      target_fasta_path = os.path.join(out_dir,"accpro_target.fasta") 
      SaveSequence(s, target_fasta_path)
      aln_filename = structanno._PrepareSS_ACCPro(p['msa'])
      _, accpro = structanno._GetSS_ACCproAnnotation(target_fasta_path,
                                                     aln_filename,
                                                     len(s))
      accpro = "".join(accpro)
      seq_to_save = CreateSequence("accpro",accpro)
      SaveSequence(seq_to_save,accpro_file,format="fasta")
      acc_files.append(accpro_file)
      CopyToCache(accpro_cache,h,accpro_file)
    except:
      accpro_failed = True

  #lets gather the final predictions
  psipred_pred = list()
  psipred_conf = list()
  accpro_pred = list()

  for i,p in enumerate(profiles):
    psipred_pred.append("".join([item for item in p["ss_pred"]]))
    psipred_conf.append("".join([str(item) for item in p["ss_conf"]]))
    if not accpro_failed:
      acc_seq = LoadSequence(acc_files[i], format = "fasta")
      accpro_pred.append(acc_seq.GetString())

  #and generate the QMEAN PSIPREDHandler/ACCPROHandler
  psipred_handler = list()
  for s,p,c in zip(sequences,psipred_pred,psipred_conf):
    data = dict()
    data["seq"] = s.GetString()
    data["ss"] = p
    data["conf"] = c
    pp = PSIPREDHandler(data)
    psipred_handler.append(pp)

  accpro_handler = None
  if not accpro_failed:
    accpro_handler = list()
    for s,a in zip(sequences,accpro_pred):
      data = dict()
      data["seq"] = s.GetString()
      data["acc"] = a
      acc = ACCPROHandler(data)
      accpro_handler.append(acc)

  return (psipred_handler, accpro_handler)


def GetDistanceConstraints(sequences, hhblits_cache):
  return None

def GetMemParam(model,working_dir):

  membrane_finder_input = model.Select('ele!=H')
  surf = msms.CalculateSurface(membrane_finder_input,radius=1.4,msms_exe=config.MSMS_BIN)[0]
  naccess.CalculateSurfaceArea(membrane_finder_input, scratch_dir = working_dir, naccess_exe=config.NACCESS_BIN)
  asa = list()
  for a in membrane_finder_input.atoms:
    if a.HasProp('asaAtom'):
      asa.append(a.GetFloatProp('asaAtom'))
    else:
      asa.append(0.0)
  mem_param = FindMembrane(mol.CreateEntityFromView(membrane_finder_input,False), surf, asa)
  return mem_param


USAGE = "usage: sm run_qmean.py <INPUT_DIR> <OUTPUT_DIR> <CACHE_DIR> <TMP_DIR>"

if len(sys.argv) < 5:
  print USAGE
  sys.exit(-1)

in_dir = sys.argv[1]
out_dir = sys.argv[2]
cache_dir = sys.argv[3]
tmp_dir = sys.argv[4]

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

#lets extract the names of the models and the sequences
model_files = list()
sequences = list()
for item in project_data["models"]:
  model_files.append(str(item["modelid"])+".pdb")
  seq_list = CreateSequenceList()
  for s in item["seqres"]:
    seq_list.AddSequence(ost.seq.CreateSequence(str(s["name"]),str(s["sequence"])))
  sequences.append(seq_list)

#iterate over all models and do the quality assessment
for i,f in enumerate(model_files):

  do_local_assessment = True
  do_global_assessment = True
  out_path = os.path.join(out_dir,f.split('.')[0])
  model = LoadPDB(os.path.join(in_dir,f)).Select("peptide=true")

  #Let's get all stuff we predict only with the sequence
  psipred_handler, accpro_handler = GetSequenceFeatures(sequences[i], 
                                                        hhblits_cache, 
                                                        accpro_cache, 
                                                        out_path)


  dc = None
  if use_qmeandisco:
    dc = GetDistanceConstraints(sequences,hhblits_cache)

  membrane_representation = None
  if use_qmeanbrane:
    mem_p = GetMemParam(model,tmp_dir) 
    membrane_representation = CreateMembraneRepresentation(model, mem_param = mem_p)
    mqa_result_membrane.AssessModelQuality(model,
                                           mem_param = mem_p,
                                           output_dir = out_path,
                                           psipred = psipred_handler,
                                           accpro = accpro_handler,
                                           dssp_path = config.DSSP_BIN)
    do_local_assessment = False
    

  mqa_result.AssessModelQuality(model,
                                local_scores = do_local_assessment,
                                global_scores = do_global_assessment, 
                                output_dir = out_path, 
                                psipred = psipred_handler,
                                accpro = accpro_handler,
                                dssp_path = config.DSSP_BIN)


  SavePDB(model,os.path.join(out_path,"model.pdb"))
  if membrane_representation != None:
    SavePDB(membrane_representation, os.path.join(out_path,"dummy_mem.pdb"))


#let's change the status to completed
status_file = open(status_path,'w')
msg = status_jobid + " " + "COMPLETED"
status_file.write(msg)
status_file.close()

