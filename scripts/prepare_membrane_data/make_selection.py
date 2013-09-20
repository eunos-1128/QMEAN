import re
import os
from opm_parser import *
from ost.bindings import clustalw
from ost import seq
from ost.bindings import kclust
from ost.bindings import dssp
from ost.table import *

#This script selects all structures, that fulfill a certain resolution
#threshold and have a minimal amount of membrane associated residues.
#It directly clusters them, calculates an according weight and stores the stuff 
#with additional information in a ost table.
#It requires as input: -The directory, where the opm structures are 
#                       (download them at: http://opm.phar.umich.edu/)
#                      -The opm file defining the transmembrane regions
#                       (download at: http://opm.phar.umich.edu/subunits.php)
#                      -A file containing resolution information in format: "id  ;  reso"
#                       on BC2: /import/bc2/data/databases/PDB/derived_data/index/resolu.idx
#                      -The name of the output table
#                      -resolution threshold
#                      -residue cutoff
#
#all parameters have to be controlled below


opm_dir='OPMStructures/'
opm_data=open('OPMTMSegments.txt').readlines()
resolu_data=open('resolu.idx').readlines()
out_tab_name='chain_weights.csv'

resolu_cutoff=2.5
residue_cutoff=30

if len(out_tab_name.split('.'))==1:
  table_format='csv'
else:
  table_format=out_tab_name.split('.')[-1]

out_tab=Table()
out_tab.AddCol('pdb_id','s')
out_tab.AddCol('chain_name','s')
out_tab.AddCol('weight','f')
out_tab.AddCol('structure_type','s')
out_tab.AddCol('membrane_residues','i')
out_tab.AddCol('resolution','f')
out_tab.AddCol('cluster','i')
out_tab.AddCol('full_query','s')

#get all ost selection queries from the provided opm file

queries=ParseOPM(opm_data)
sequences=seq.CreateSequenceList()
structure_type_info=dict()
residue_number_info=dict()
resolution_info=dict()

fail_ids=list()

for id, query in queries.iteritems():

  fulfills_resolution=False
  fulfills_residues=False

  #check wether structure fulfills resolution threshold
  for l in resolu_data:
    if l[:4].upper()==id.upper():
      res=float(l.split()[2])
      if res<=resolu_cutoff and res>0:
        fulfills_resolution=True
        resolution_info[id]=res
      break

  #check wether enough residues are in the membrane
  number_of_residues=0
  regex=r'(?P<start>\d+):(?P<end>\d+)'

  for segment in re.findall(regex,query):
    number_of_residues+=(int(segment[1])-int(segment[0])+1)

  if number_of_residues > residue_cutoff:
    fulfills_residues=True

  if fulfills_resolution==False or fulfills_residues==False:
    continue

  try:
    prot=io.LoadPDB(opm_dir+id+'.pdb').Select('peptide=true')
  except:
    print 'could not load '+id
    fail_ids.append(id)
    continue

  dssp.AssignDSSP(prot)

  membrane_view=prot.Select(query)

  alpha_residues=0
  beta_residues=0

  #get structure type

  for r in membrane_view.residues:
    if r.GetSecStructure().IsHelical():
      alpha_residues+=1
    if r.GetSecStructure().IsExtended():
      beta_residues+=1

  structure_type=''

  if alpha_residues>beta_residues:
    structure_type='alpha'
  else:
    structure_type='beta'

  structure_type_info[id]=structure_type
      
  #Add chains, that are membrane associated. Note, that the sequence
  #of the FULL chain (with soluble part) is used for clustering.
  #In case of homooligomers, only one chain is added.

  used_chain_sequences=list()

  for c in membrane_view.chains:
    s=''
    for i,r in enumerate(c.handle.Select('peptide=true').residues):
      s+=r.one_letter_code

    new_seq=seq.CreateSequence(id+'_'+c.name,s)
    already_present=False

    for sequence in used_chain_sequences:
      try:
        aln=clustalw.ClustalW(sequence,new_seq)
      except:
        already_present=True
        continue
      if seq.alg.SequenceIdentity(aln)>90:
        already_present=True

    if already_present:
      continue
    else:
      used_chain_sequences.append(new_seq)
      residue_number_info[id+'_'+c.name]=len(c.residues)

  for s in used_chain_sequences:
    sequences.AddSequence(s)

kclustout=kclust.kClust(sequences)

for i,c in enumerate(kclustout):
  sequences=c.sequences
  for s in sequences:
    row=dict()
    row['pdb_id']=s.GetName()[:4]
    row['chain_name']=s.GetName()[-1]
    row['weight']=1.0/len(sequences)
    row['structure_type']=structure_type_info[s.GetName()[:4]]
    row['membrane_residues']=residue_number_info[s.GetName()]
    row['resolution']=resolution_info[s.GetName()[:4]]
    row['cluster']=i
    row['full_query']=queries[s.GetName()[:4]]
    out_tab.AddRow(row)

out_tab.Save(out_tab_name, format=table_format)
print 'could not load following structures:'
for f_id in fail_ids:
  print f_id
