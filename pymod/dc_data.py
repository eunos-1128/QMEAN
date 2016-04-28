import qmean
import numpy as np
import scipy as sp
from scipy import cluster,spatial
try:
  import sm
  from sm.pipeline import*
  from sm.smtl import *
  from sm.core import hhbase
except:
  print "Could not import SWISS-MODEL"  
from operator import itemgetter
import os, shutil, md5
from ost.seq.alg import BLOSUM62, SequenceSimilarity


class TemplateInformation:
  def __init__(self, seqres, working_dir = '.'):
    if not os.path.exists(working_dir):
      os.makedirs(working_dir)
      print "Created directory: %s" %working_dir
    self.working_dir = working_dir
    self.seqres = seqres

  def FindTemplates(self):  
    # creation of Swissmodel project and template search
    md5_hash = md5.new(self.seqres).hexdigest()
    prj_name = md5_hash
    sm_prj_path = self.working_dir+'/'+prj_name
    cmd = 'sm init '+sm_prj_path+' '+self.seqres

    os.system(cmd)
    cmd  = 'sm tplsearch '+sm_prj_path+'.sm --hhblits --no-ostate-prediction'
    os.system(cmd)
    prj = Project(sm_prj_path+'.sm')
    tpls = prj.GetTemplates()

    # align template sequences with seqres and find SMTL filename of templates
    aln_dict = dict()      

    if len(tpls) > 65500:
      print "Too many templates. Only first 65500 templates will be considered."
      tpls = tpls[:65500];

    for tpl in tpls:
      pdb_id = tpl.pdb_id
      ass_id = tpl.assembly_id
      chain_name = tpl.chain_name

      # prune alignment of target sequence with template sequence   
      bio_unit = prj.smtl.Get(pdb_id,ass_id)
      new_aln = HomologyModel.PruneAtomSeqAlignment(bio_unit.GetChainByName(chain_name).ToAtomSeqAlignment(tpl.alignment),2)        
      fn = prj.GetModelRepos().FilenameForModel(pdb_id,ass_id,chain_name)
      aln_dict[fn] = new_aln 

    # remove Swissmodel project
    shutil.rmtree(sm_prj_path+'.sm')     

    # multiple sequence alignment  
    alns = ost.seq.AlignmentList()
    tpl_list = list()
    for k,aln in sorted(aln_dict.items()):
      alns.append(aln)
      tpl_list.append(k)  
    ref_seq = seq.CreateSequence("ref", self.seqres)
    msaln = ost.seq.alg.MergePairwiseAlignments(alns,ref_seq)

    return msaln, tpl_list


  def ClusterTemplates(self,msaln):
    weight_mat= BLOSUM62
    num_tpls = msaln.GetCount()-1
    sim_mat= np.empty([num_tpls,num_tpls])
    
    # no templates  
    if num_tpls == 1: 
      return []
    # only one template. A list with only one cluster is returned  
    if num_tpls == 2: 
      return [[0]]  

    for i in range(1,msaln.GetCount()):
      for j in range(i,msaln.GetCount()):
        actual_sim = SequenceSimilarity(msaln,weight_mat,False,i,j)
        sim_mat[i-1,j-1] = actual_sim
        sim_mat[j-1,i-1] = actual_sim

    # cluster templates by sequence similarity
    condensed = sp.spatial.distance.pdist(sim_mat)
    h_cluster = sp.cluster.hierarchy.linkage(condensed) 
    f_cluster = sp.cluster.hierarchy.fcluster(h_cluster,2,criterion= 'distance')
    assignment = sp.array(sorted(zip(list(f_cluster),range(num_tpls)),key= itemgetter(0)))

    split_list = list()
    for i,ele in enumerate(assignment):
      if ele[0]!=assignment[i-1][0]:
        split_list.append(i)

    clusters = sp.split(assignment,split_list[1:])
    cluster_list = list()
    for cl in clusters:
      cluster_list.append(cl[:,1].tolist())

    return cluster_list
