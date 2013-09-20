from qmean import *
from ost.table import *
from ost.bindings import dssp
import os
import traceback
import ost


def BuildSolubleInteraction(lower_cutoff, upper_cutoff, dist_bins, seq_sep, pdb_dir, id_information, ref_state='classic',s=0.02):
  """
  returns a dictionary with keys ['stat','pot'], that contains statistic, potentialcontainer respectively built
  with the given parameters. The statistics/potentials in the conatiner have the names ['helix','extended','coil']
  lower_cutoff: minimum distance in interaction 
  upper_cutoff: maximum distance in interaction
  dist_bins: number of distance bins 
  seq_sep: sequence separation
  pdb_dir: directory in which the training pdbs are stored
  id_information: a dictionary with the filename of the structures as key and a chain, from which you want to extract
                  the potentials as value. The filenames must be present in the given pdb_dir
  ref_state: the reference state, that will be used for the creation of the potentials out of the statistics 
  s: sigma parameter important for handling sparse data. See "Calculation of Conformational Ensembles from Potentials of Mean Force" M.J. Sippl 1990

  """

  ret_dict=dict()

  stat=StatisticContainer()
  pot=PotentialContainer()

  stat['helix']=InteractionStatistic(lower_cutoff, upper_cutoff, dist_bins, seq_sep)
  stat['extended']=InteractionStatistic(lower_cutoff, upper_cutoff, dist_bins, seq_sep)  
  stat['coil']=InteractionStatistic(lower_cutoff, upper_cutoff, dist_bins, seq_sep)

  num_alpha=0
  num_beta=0
  num_coil=0

  for k,v in id_information.iteritems():
    try:

      print 'extracting statistics from ', k
      ent=ost.io.LoadPDB(os.path.join(pdb_dir,k+'.pdb')).Select('peptide=true and occ=1')

      dssp.AssignDSSP(ent)

      alpha_selection=ent.Select('chain='+v+' and rtype=helix')
      beta_selection=ent.Select('chain='+v+' and rtype=ext')
      coil_selection=ent.Select('chain='+v+' and rtype=coil')

      stat['helix'].Extract(alpha_selection,ent)
      stat['extended'].Extract(beta_selection,ent)
      stat['coil'].Extract(coil_selection,ent)

      num_alpha+=len(alpha_selection.residues)
      num_beta+=len(beta_selection.residues)
      num_coil+=len(coil_selection.residues)

    except:
      print traceback.print_exc()
      pass

  pot['helix']=InteractionPotential.Create(stat['helix'], sigma=s, reference_state=ref_state)
  pot['extended']=InteractionPotential.Create(stat['extended'], sigma=s, reference_state=ref_state)
  pot['coil']=InteractionPotential.Create(stat['coil'], sigma=s, reference_state=ref_state)

  ret_dict['stat']=stat
  ret_dict['pot']=pot
  ret_dict['statistics']={'num_alpha':num_alpha,'num_beta':num_beta,'num_coil':num_coil}
  return ret_dict


def BuildSolubleCBeta(lower_cutoff, upper_cutoff, dist_bins, seq_sep, pdb_dir, id_information, ref_state='classic',s=0.02):
  """
  returns a dictionary with keys ['stat','pot'], that contains statistic, potentialcontainer respectively built
  with the given parameters. The statistics/potentials in the conatiner have the names ['helix','extended','coil']
  lower_cutoff: minimum distance in cbeta 
  upper_cutoff: maximum distance in cbeta
  dist_bins: number of distance bins 
  seq_sep: sequence separation
  pdb_dir: directory in which the training pdbs are stored
  id_information: a dictionary with the filename of the structures as key and a chain, from which you want to extract
                  the potentials as value. The filenames must be present in the given pdb_dir
  ref_state: the reference state, that will be used for the creation of the potentials out of the statistics
  s: sigma parameter important for handling sparse data. See "Calculation of Conformational Ensembles from Potentials of Mean Force" M.J. Sippl 1990

  """

  ret_dict=dict()

  stat=StatisticContainer()
  pot=PotentialContainer()

  stat['helix']=CBetaStatistic(lower_cutoff, upper_cutoff, dist_bins, seq_sep)
  stat['extended']=CBetaStatistic(lower_cutoff, upper_cutoff, dist_bins, seq_sep)  
  stat['coil']=CBetaStatistic(lower_cutoff, upper_cutoff, dist_bins, seq_sep)

  num_alpha=0
  num_beta=0
  num_coil=0

  for k,v in id_information.iteritems():
    try:

      print 'extracting statistics from ', k

      ent=ost.io.LoadPDB(os.path.join(pdb_dir,k+'.pdb')).Select('peptide=true and occ=1')

      dssp.AssignDSSP(ent)

      alpha_selection=ent.Select('chain='+v+' and rtype=helix')
      beta_selection=ent.Select('chain='+v+' and rtype=ext')
      coil_selection=ent.Select('chain='+v+' and rtype=coil')

      stat['helix'].Extract(alpha_selection,ent)
      stat['extended'].Extract(beta_selection,ent)
      stat['coil'].Extract(coil_selection,ent)

      num_alpha+=len(alpha_selection.residues)
      num_beta+=len(beta_selection.residues)
      num_coil+=len(coil_selection.residues)

    except:
      pass

  pot['helix']=CBetaPotential.Create(stat['helix'], sigma=s, reference_state=ref_state)
  pot['extended']=CBetaPotential.Create(stat['extended'], sigma=s, reference_state=ref_state)
  pot['coil']=CBetaPotential.Create(stat['coil'], sigma=s, reference_state=ref_state)

  ret_dict['stat']=stat
  ret_dict['pot']=pot
  ret_dict['statistics']={'num_alpha':num_alpha,'num_beta':num_beta,'num_coil':num_coil}
  return ret_dict


def BuildSolublePacking(cutoff, max_count, bin_size, pdb_dir, id_information, ref_state='classic',s=0.02):
  """
  returns a dictionary with keys ['stat','pot'], that contains statistic, potentialcontainer respectively built
  with the given parameters. The statistics/potentials in the conatiner have the names ['helix','extended','coil']
  lower_cutoff: minimum distance in interaction 
  max_count: upper bound for counting atoms in the neighbourhood
  bin_size: size of bin in atom counts per bin  
  pdb_dir: directory in which the training pdbs are stored
  id_information: a dictionary with the filename of the structures as key and a chain, from which you want to extract
                  the potentials as value. The filenames must be present in the given pdb_dir
  ref_state: the reference state, that will be used for the creation of the potentials out of the statistics
  s: sigma parameter important for handling sparse data. See "Calculation of Conformational Ensembles from Potentials of Mean Force" M.J. Sippl 1990

  """

  ret_dict=dict()

  stat=StatisticContainer()
  pot=PotentialContainer()

  stat['helix']=PackingStatistic(cutoff, max_count, bin_size)
  stat['extended']=PackingStatistic(cutoff, max_count, bin_size)  
  stat['coil']=PackingStatistic(cutoff, max_count, bin_size)

  num_alpha=0
  num_beta=0
  num_coil=0

  for k,v in id_information.iteritems():
    try:

      print 'extracting statistics from ', k

      ent=ost.io.LoadPDB(os.path.join(pdb_dir,k+'.pdb')).Select('peptide=true')

      dssp.AssignDSSP(ent)

      alpha_selection=ent.Select('chain='+v+' and rtype=helix')
      beta_selection=ent.Select('chain='+v+' and rtype=ext')
      coil_selection=ent.Select('chain='+v+' and rtype=coil')

      stat['helix'].Extract(alpha_selection,ent)
      stat['extended'].Extract(beta_selection,ent)
      stat['coil'].Extract(coil_selection,ent)

      num_alpha+=len(alpha_selection.residues)
      num_beta+=len(beta_selection.residues)
      num_coil+=len(coil_selection.residues)

    except:
      pass

  pot['helix']=PackingPotential.Create(stat['helix'], sigma=s, reference_state=ref_state)
  pot['extended']=PackingPotential.Create(stat['extended'], sigma=s, reference_state=ref_state)
  pot['coil']=PackingPotential.Create(stat['coil'], sigma=s, reference_state=ref_state)

  ret_dict['stat']=stat
  ret_dict['pot']=pot
  ret_dict['statistics']={'num_alpha':num_alpha,'num_beta':num_beta,'num_coil':num_coil}
  return ret_dict



def BuildSolubleTorsion(group_identifier, angle_bins, pdb_dir, id_information, ref_state='classic',s=0.02):
  """
  returns a dictionary with keys ['stat','pot'], that contains statistic, potentialcontainer respectively built
  with the given parameters. The statistic/potential in the conatiner have the names 'overall'
  group_identifier: a list of strings, that describes the different groups of the torsion potential.
                    e.g. 'ALA,PHE-HIS-PRO'
  angle_bins: number of bins for the corresponding angles. The parameter must be a list with 6 entries
              =>[phi_prev, psi_prev, phi_central, psi_central, phi_next, psi_next] 
  pdb_dir: directory in which the training pdbs are stored
  id_information: a dictionary with the filename of the structures as key and a chain, from which you want to extract
                  the potentials as value. The filenames must be present in the given pdb_dir
  ref_state: the reference state, that will be used for the creation of the potentials out of the statistics
  s: sigma parameter important for handling sparse data. See "Calculation of Conformational Ensembles from Potentials of Mean Force" M.J. Sippl 1990

  """

  ret_dict=dict()

  stat=StatisticContainer()
  pot=PotentialContainer()

  stat['overall']=TorsionStatistic(group_identifier, angle_bins)

  for k,v in id_information.iteritems():
    try:

      print 'extracting statistics from ', k

      ent=ost.io.LoadPDB(os.path.join(pdb_dir,k+'.pdb')).Select('peptide=true')

      stat['overall'].Extract(ent)


    except:
      pass

  pot['overall']=TorsionPotential.Create(stat['overall'], sigma=s, reference_state=ref_state)

  ret_dict['stat']=stat
  ret_dict['pot']=pot

  return ret_dict



def BuildSolubleTorsionSeparateSecondaryStructures(group_identifier, angle_bins, pdb_dir, id_information, ref_state='classic',s=0.02):
  """
  returns a dictionary with keys ['stat','pot'], that contains statistic, potentialcontainer respectively built
  with the given parameters. The statistics/potentials in the conatiner have the names ['helix','extended','coil']
  group_identifier: a list of strings, that describes the different groups of the torsion potential.
                    e.g. 'ALA,PHE-HIS-PRO'
  angle_bins: number of bins for the corresponding angles. The parameter must be a list with 6 entries
              =>[phi_prev, psi_prev, phi_central, psi_central, phi_next, psi_next] 
  pdb_dir: directory in which the training pdbs are stored
  id_information: a dictionary with the filename of the structures as key and a chain, from which you want to extract
                  the potentials as value. The filenames must be present in the given pdb_dir
  ref_state: the reference state, that will be used for the creation of the potentials out of the statistics
  s: sigma parameter important for handling sparse data. See "Calculation of Conformational Ensembles from Potentials of Mean Force" M.J. Sippl 1990

  """

  ret_dict=dict()

  stat=StatisticContainer()
  pot=PotentialContainer()

  num_alpha=0
  num_beta=0
  num_coil=0

  stat['helix']=TorsionStatistic(group_identifier, angle_bins)
  stat['extended']=TorsionStatistic(group_identifier, angle_bins)
  stat['coil']=TorsionStatistic(group_identifier, angle_bins)

  for k,v in id_information.iteritems():
    try:

      print 'extracting statistics from ', k

      ent=ost.io.LoadPDB(os.path.join(pdb_dir,k+'.pdb')).Select('peptide=true')

      dssp.AssignDSSP(ent)

      alpha_selection=ent.Select('chain='+v+' and rtype=helix')
      beta_selection=ent.Select('chain='+v+' and rtype=ext')
      coil_selection=ent.Select('chain='+v+' and rtype=coil')

      stat['helix'].Extract(alpha_selection)
      stat['extended'].Extract(beta_selection)
      stat['coil'].Extract(coil_selection)

      num_alpha+=len(alpha_selection.residues)
      num_beta+=len(beta_selection.residues)
      num_coil+=len(coil_selection.residues)

    except:
      pass

  pot['helix']=TorsionPotential.Create(stat['helix'], sigma=s, reference_state=ref_state)
  pot['extended']=TorsionPotential.Create(stat['extended'], sigma=s, reference_state=ref_state)
  pot['coil']=TorsionPotential.Create(stat['coil'], sigma=s, reference_state=ref_state)

  ret_dict['stat']=stat
  ret_dict['pot']=pot
  ret_dict['statistics']={'num_alpha':num_alpha,'num_beta':num_beta,'num_coil':num_coil}
  return ret_dict



def BuildSolubleReduced(lower_cutoff, upper_cutoff, angle_bins, dist_bins, seq_sep, pdb_dir, id_information, ref_state='classic',s=0.02):
  """
  returns a dictionary with keys ['stat','pot'], that contains statistic, potentialcontainer respectively built
  with the given parameters. The statistics/potentials in the conatiner have the names ['helix','extended','coil']
  lower_cutoff: minimum CA-distance in reduced 
  upper_cutoff: maximum CA-distance in reduced
  dist_bins: number of distance bins
  angle_bins: number of angular bins 
  seq_sep: sequence separation
  pdb_dir: directory in which the training pdbs are stored
  id_information: a dictionary with the filename of the structures as key and a chain, from which you want to extract
                  the potentials as value. The filenames must be present in the given pdb_dir
  ref_state: the reference state, that will be used for the creation of the potentials out of the statistics
  s: sigma parameter important for handling sparse data. See "Calculation of Conformational Ensembles from Potentials of Mean Force" M.J. Sippl 1990

  """

  ret_dict=dict()

  stat=StatisticContainer()
  pot=PotentialContainer()

  stat['helix']=ReducedStatistic(lower_cutoff, upper_cutoff, angle_bins, dist_bins, seq_sep)
  stat['extended']=ReducedStatistic(lower_cutoff, upper_cutoff, angle_bins, dist_bins, seq_sep)  
  stat['coil']=ReducedStatistic(lower_cutoff, upper_cutoff, angle_bins, dist_bins, seq_sep)

  num_alpha=0
  num_beta=0
  num_coil=0

  for k,v in id_information.iteritems():
    try:

      print 'extracting statistics from ', k

      ent=ost.io.LoadPDB(os.path.join(pdb_dir,k+'.pdb')).Select('peptide=true')

      dssp.AssignDSSP(ent)

      alpha_selection=ent.Select('chain='+v+' and rtype=helix')
      beta_selection=ent.Select('chain='+v+' and rtype=ext')
      coil_selection=ent.Select('chain='+v+' and rtype=coil')

      stat['helix'].Extract(alpha_selection,ent)
      stat['extended'].Extract(beta_selection,ent)
      stat['coil'].Extract(coil_selection,ent)

      num_alpha+=len(alpha_selection.residues)
      num_beta+=len(beta_selection.residues)
      num_coil+=len(coil_selection.residues)

    except:
      pass

  pot['helix']=ReducedPotential.Create(stat['helix'], sigma=s, reference_state=ref_state)
  pot['extended']=ReducedPotential.Create(stat['extended'], sigma=s, reference_state=ref_state)
  pot['coil']=ReducedPotential.Create(stat['coil'], sigma=s, reference_state=ref_state)

  ret_dict['stat']=stat
  ret_dict['pot']=pot
  ret_dict['statistics']={'num_alpha':num_alpha,'num_beta':num_beta,'num_coil':num_coil}
  return ret_dict






def BuildMembraneInteraction(lower_cutoff, upper_cutoff, dist_bins, seq_sep, pdb_dir, info_table, ref_state='classic',s=0.02):
  """
  returns a dictionary with keys ['stat','pot'], that contains statistic, potentialcontainer respectively built
  with the given parameters. The statistics/potentials in the conatiner have the names ['helix','extended','coil']
  lower_cutoff: minimum distance in interaction 
  upper_cutoff: maximum distance in interaction
  dist_bins: number of distance bins 
  seq_sep: sequence separation
  pdb_dir: directory in which the training pdbs are stored
  info_table: output of the make_selection.py script in qmean/scripts/prepare_membrane_data. An ost table containing at least following columns:
              ['pdb_id', 'chain_name', 'weight', 'structure_type', full_query']
              pdb_id: as it says...
              chain_name: as it says...
              weight: value in range [0.0, 1.0]
              structure_type: structure type of transmembrane part can either be 'alpha' or 'beta'
              full_query: ost query of membrane part of the full structure


  ref_state: the reference state, that will be used for the creation of the potentials out of the statistics
  s: sigma parameter important for handling sparse data. See "Calculation of Conformational Ensembles from Potentials of Mean Force" M.J. Sippl 1990

  """

  ret_dict=dict()

  stat=StatisticContainer()
  pot=PotentialContainer()

  num_alpha=0
  num_beta=0

  stat['helix']=InteractionStatistic(lower_cutoff, upper_cutoff, dist_bins, seq_sep)
  stat['extended']=InteractionStatistic(lower_cutoff, upper_cutoff, dist_bins, seq_sep)  

  id_idx=info_table.GetColIndex('pdb_id')
  chain_idx=info_table.GetColIndex('chain_name')
  weight_idx=info_table.GetColIndex('weight')
  structure_type_idx=info_table.GetColIndex('structure_type')
  query_idx=info_table.GetColIndex('full_query')

  for r in info_table.rows:
    try:

      print 'extracting statistics from ', r[id_idx]

      ent=ost.io.LoadPDB(os.path.join(pdb_dir,r[id_idx]+'.pdb')).Select('peptide=true')
      membrane_part=ent.Select(r[query_idx])
      chain_selection_membrane=membrane_part.Select('cname='+r[chain_idx])
      
      if r[structure_type_idx]=='alpha':
        stat['helix'].Extract(chain_selection_membrane, membrane_part, weight=r[weight_idx])
        num_alpha+=len(chain_selection_membrane.residues)
      elif r[structure_type_idx]=='beta':
        stat['extended'].Extract(chain_selection_membrane, membrane_part, weight=r[weight_idx])
        num_beta+=len(chain_selection_membrane.residues)

    except:
      print traceback.print_exc()

  pot['helix']=InteractionPotential.Create(stat['helix'], sigma=s, reference_state=ref_state)
  pot['extended']=InteractionPotential.Create(stat['extended'], sigma=s, reference_state=ref_state)

  ret_dict['stat']=stat
  ret_dict['pot']=pot
  ret_dict['statistics']={'num_alpha':num_alpha,'num_beta':num_beta}

  return ret_dict


def BuildMembraneCBeta(lower_cutoff, upper_cutoff, dist_bins, seq_sep, pdb_dir, info_table, ref_state='classic',s=0.02):
  """
  returns a dictionary with keys ['stat','pot'], that contains statistic, potentialcontainer respectively built
  with the given parameters. The statistics/potentials in the conatiner have the names ['helix','extended','coil']
  lower_cutoff: minimum distance in cbeta 
  upper_cutoff: maximum distance in cbeta
  dist_bins: number of distance bins 
  seq_sep: sequence separation
  pdb_dir: directory in which the training pdbs are stored
  info_table: output of the make_selection.py script in qmean/scripts/prepare_membrane_data. An ost table containing at least following columns:
              ['pdb_id', 'chain_name', 'weight', 'structure_type', full_query']
              pdb_id: as it says...
              chain_name: as it says...
              weight: value in range [0.0, 1.0]
              structure_type: structure type of transmembrane part can either be 'alpha' or 'beta'
              full_query: ost query of membrane part of the full structure

  ref_state: the reference state, that will be used for the creation of the potentials out of the statistics
  s: sigma parameter important for handling sparse data. See "Calculation of Conformational Ensembles from Potentials of Mean Force" M.J. Sippl 1990

  """

  ret_dict=dict()

  stat=StatisticContainer()
  pot=PotentialContainer()

  num_alpha=0
  num_beta=0

  stat['helix']=CBetaStatistic(lower_cutoff, upper_cutoff, dist_bins, seq_sep)
  stat['extended']=CBetaStatistic(lower_cutoff, upper_cutoff, dist_bins, seq_sep)

  id_idx=info_table.GetColIndex('pdb_id')
  chain_idx=info_table.GetColIndex('chain_name')
  weight_idx=info_table.GetColIndex('weight')
  structure_type_idx=info_table.GetColIndex('structure_type')
  query_idx=info_table.GetColIndex('full_query')

  for r in info_table.rows:
    try:

      print 'extracting statistics from ', r[id_idx]

      ent=ost.io.LoadPDB(os.path.join(pdb_dir,r[id_idx]+'.pdb')).Select('peptide=true')
      membrane_part=ent.Select(r[query_idx])
      chain_selection_membrane=membrane_part.Select('cname='+r[chain_idx])
      
      if r[structure_type_idx]=='alpha':
        stat['helix'].Extract(chain_selection_membrane, membrane_part, weight=r[weight_idx])
        num_alpha+=len(chain_selection_membrane.residues)
      elif r[structure_type_idx]=='beta':
        stat['extended'].Extract(chain_selection_membrane, membrane_part, weight=r[weight_idx])
        num_beta+=len(chain_selection_membrane.residues)

    except:
      pass

  pot['helix']=CBetaPotential.Create(stat['helix'], sigma=s, reference_state=ref_state)
  pot['extended']=CBetaPotential.Create(stat['extended'], sigma=s, reference_state=ref_state)


  ret_dict['stat']=stat
  ret_dict['pot']=pot
  ret_dict['statistics']={'num_alpha':num_alpha,'num_beta':num_beta}
  return ret_dict


def BuildMembranePacking(cutoff, max_count, bin_size, pdb_dir, info_table, ref_state='classic',s=0.02):
  """
  returns a dictionary with keys ['stat','pot'], that contains statistic, potentialcontainer respectively built
  with the given parameters. The statistics/potentials in the conatiner have the names ['helix','extended','coil']
  lower_cutoff: minimum distance in interaction 
  max_count: upper bound for counting atoms in the neighbourhood
  bin_size: size of bin in atom counts per bin  
  pdb_dir: directory in which the training pdbs are stored
  info_table: output of the make_selection.py script in qmean/scripts/prepare_membrane_data. An ost table containing at least following columns:
              ['pdb_id', 'chain_name', 'weight', 'structure_type', full_query']
              pdb_id: as it says...
              chain_name: as it says...
              weight: value in range [0.0, 1.0]
              structure_type: structure type of transmembrane part can either be 'alpha' or 'beta'
              full_query: ost query of membrane part of the full structure

  ref_state: the reference state, that will be used for the creation of the potentials out of the statistics
  s: sigma parameter important for handling sparse data. See "Calculation of Conformational Ensembles from Potentials of Mean Force" M.J. Sippl 1990

  """

  ret_dict=dict()

  stat=StatisticContainer()
  pot=PotentialContainer()

  num_alpha=0
  num_beta=0

  stat['helix']=PackingStatistic(cutoff, max_count, bin_size)
  stat['extended']=PackingStatistic(cutoff, max_count, bin_size)  


  id_idx=info_table.GetColIndex('pdb_id')
  chain_idx=info_table.GetColIndex('chain_name')
  weight_idx=info_table.GetColIndex('weight')
  structure_type_idx=info_table.GetColIndex('structure_type')
  query_idx=info_table.GetColIndex('full_query')

  for r in info_table.rows:
    try:

      print 'extracting statistics from ', r[id_idx]

      ent=ost.io.LoadPDB(os.path.join(pdb_dir,r[id_idx]+'.pdb')).Select('peptide=true')
      membrane_part=ent.Select(r[query_idx])
      chain_selection_membrane=membrane_part.Select('cname='+r[chain_idx])
      
      if r[structure_type_idx]=='alpha':
        stat['helix'].Extract(chain_selection_membrane, ent, weight=r[weight_idx])
        num_alpha+=len(chain_selection_membrane.residues)
      elif r[structure_type_idx]=='beta':
        stat['extended'].Extract(chain_selection_membrane, ent, weight=r[weight_idx])
        num_beta+=len(chain_selection_membrane.residues)

    except:
      pass

  pot['helix']=PackingPotential.Create(stat['helix'], sigma=s, reference_state=ref_state)
  pot['extended']=PackingPotential.Create(stat['extended'], sigma=s, reference_state=ref_state)

  ret_dict['stat']=stat
  ret_dict['pot']=pot
  ret_dict['statistics']={'num_alpha':num_alpha,'num_beta':num_beta}
  return ret_dict


def BuildMembraneTorsion(group_identifier, angle_bins, pdb_dir, info_table, ref_state='classic',s=0.02):
  """
  returns a dictionary with keys ['stat','pot'], that contains statistic, potentialcontainer respectively built
  with the given parameters. The statistic/potential in the conatiner have the names 'overall'
  group_identifier: a list of strings, that describes the different groups of the torsion potential.
                    e.g. 'ALA,PHE-HIS-PRO'
  angle_bins: number of bins for the corresponding angles. The parameter must be a list with 6 entries
              =>[phi_prev, psi_prev, phi_central, psi_central, phi_next, psi_next] 
  pdb_dir: directory in which the training pdbs are stored
  info_table: output of the make_selection.py script in qmean/scripts/prepare_membrane_data. An ost table containing at least following columns:
              ['pdb_id', 'chain_name', 'weight', 'structure_type', full_query']
              pdb_id: as it says...
              chain_name: as it says...
              weight: value in range [0.0, 1.0]
              structure_type: structure type of transmembrane part can either be 'alpha' or 'beta'
              full_query: ost query of membrane part of the full structure

  ref_state: the reference state, that will be used for the creation of the potentials out of the statistics
  s: sigma parameter important for handling sparse data. See "Calculation of Conformational Ensembles from Potentials of Mean Force" M.J. Sippl 1990

  """

  ret_dict=dict()

  stat=StatisticContainer()
  pot=PotentialContainer()

  num_alpha=0
  num_beta=0

  stat['helix']=TorsionStatistic(group_identifier, angle_bins)
  stat['extended']=TorsionStatistic(group_identifier, angle_bins)

  id_idx=info_table.GetColIndex('pdb_id')
  chain_idx=info_table.GetColIndex('chain_name')
  weight_idx=info_table.GetColIndex('weight')
  structure_type_idx=info_table.GetColIndex('structure_type')
  query_idx=info_table.GetColIndex('full_query')

  for r in info_table.rows:
    try:

      print 'extracting statistics from ', r[id_idx]

      ent=ost.io.LoadPDB(os.path.join(pdb_dir,r[id_idx]+'.pdb')).Select('peptide=true')
      membrane_part=ent.Select(r[query_idx])
      chain_selection_membrane=membrane_part.Select('cname='+r[chain_idx])
      
      if r[structure_type_idx]=='alpha':
        stat['helix'].Extract(chain_selection_membrane, weight=r[weight_idx])
        num_alpha+=len(chain_selection_membrane.residues)
      elif r[structure_type_idx]=='beta':
        stat['extended'].Extract(chain_selection_membrane, weight=r[weight_idx])
        num_beta+=len(chain_selection_membrane.residues)

    except:
      pass


  pot['helix']=TorsionPotential.Create(stat['helix'], sigma=s, reference_state=ref_state)
  pot['extended']=TorsionPotential.Create(stat['extended'], sigma=s, reference_state=ref_state)

  ret_dict['stat']=stat
  ret_dict['pot']=pot
  ret_dict['statistics']={'num_alpha':num_alpha,'num_beta':num_beta}
  return ret_dict



def BuildMembraneReduced(lower_cutoff, upper_cutoff, angle_bins, dist_bins, seq_sep, pdb_dir, info_table, ref_state='classic',s=0.02):
  """
  returns a dictionary with keys ['stat','pot'], that contains statistic, potentialcontainer respectively built
  with the given parameters. The statistics/potentials in the conatiner have the names ['helix','extended','coil']
  lower_cutoff: minimum CA-distance in reduced 
  upper_cutoff: maximum CA-distance in reduced
  dist_bins: number of distance bins
  angle_bins: number of angular bins 
  seq_sep: sequence separation
  pdb_dir: directory in which the training pdbs are stored
  info_table: output of the make_selection.py script in qmean/scripts/prepare_membrane_data. An ost table containing at least following columns:
              ['pdb_id', 'chain_name', 'weight', 'structure_type', full_query']
              pdb_id: as it says...
              chain_name: as it says...
              weight: value in range [0.0, 1.0]
              structure_type: structure type of transmembrane part can either be 'alpha' or 'beta'
              full_query: ost query of membrane part of the full structure

  ref_state: the reference state, that will be used for the creation of the potentials out of the statistics
  s: sigma parameter important for handling sparse data. See "Calculation of Conformational Ensembles from Potentials of Mean Force" M.J. Sippl 1990

  """

  ret_dict=dict()

  stat=StatisticContainer()
  pot=PotentialContainer()

  num_alpha=0
  num_beta=0

  stat['helix']=ReducedStatistic(lower_cutoff, upper_cutoff, angle_bins, dist_bins, seq_sep)
  stat['extended']=ReducedStatistic(lower_cutoff, upper_cutoff, angle_bins, dist_bins, seq_sep)  

  id_idx=info_table.GetColIndex('pdb_id')
  chain_idx=info_table.GetColIndex('chain_name')
  weight_idx=info_table.GetColIndex('weight')
  structure_type_idx=info_table.GetColIndex('structure_type')
  query_idx=info_table.GetColIndex('full_query')

  for r in info_table.rows:
    try:

      print 'extracting statistics from ', r[id_idx]

      ent=ost.io.LoadPDB(os.path.join(pdb_dir,r[id_idx]+'.pdb')).Select('peptide=true')
      membrane_part=ent.Select(r[query_idx])
      chain_selection_membrane=membrane_part.Select('cname='+r[chain_idx])
      
      if r[structure_type_idx]=='alpha':
        stat['helix'].Extract(chain_selection_membrane, membrane_part, weight=r[weight_idx])
        num_alpha+=len(chain_selection_membrane.residues)
      elif r[structure_type_idx]=='beta':
        stat['extended'].Extract(chain_selection_membrane, membrane_part, weight=r[weight_idx])
        num_beta+=len(chain_selection_membrane.residues)

    except:
      pass

  pot['helix']=ReducedPotential.Create(stat['helix'], sigma=s, reference_state=ref_state)
  pot['extended']=ReducedPotential.Create(stat['extended'], sigma=s, reference_state=ref_state)


  ret_dict['stat']=stat
  ret_dict['pot']=pot
  ret_dict['statistics']={'num_alpha':num_alpha,'num_beta':num_beta}
  return ret_dict
