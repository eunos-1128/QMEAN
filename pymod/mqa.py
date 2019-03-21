# Copyright (c) 2013-2018, SIB - Swiss Institute of Bioinformatics and
# Biozentrum - University of Basel
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
# http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from qmean import *
from ost import mol
from predicted_sequence_features import AlignChainToSEQRES


class Scores:

  def __init__(self, target, environment, potential_container, smooth_std=None, 
               psipred=None, accpro=None, dc=None, norm=True, count_radius=15.0):

    self.data=dict()

    #set structural data
    if isinstance(target, mol.EntityHandle):
      self.target = target.Select('')
    else:
      self.target=target
    if isinstance(environment, mol.EntityHandle):
      self.environment = environment.Select('')
    else:
      self.environment=environment

    #set predicted sequence data
    if psipred != None:
      if isinstance(psipred,list):
        if len(psipred) != len(target.chains):
          raise ValueError("List of provided PSIPREDHandlers must have same length as chains in target!")
        self.psipred = psipred
      else:
        self.psipred = len(self.target.chains) * [psipred]
    else:
      self.psipred = None

    if accpro != None:
      if isinstance(accpro,list):
        if len(accpro) != len(target.chains):
          raise ValueError("List of provided ACCPROHandlers must have same length as chains in target!")
        self.accpro = accpro
      else:
        self.accpro = len(self.target.chains) * [accpro]
    else:
      self.accpro = None

    if dc != None:  
      if isinstance(dc,list):
        if len(dc) != len(target.chains):
          raise ValueError("List of provided DCData must have same length as chains in target!")
        self.dc = dc
      else:
        self.dc = len(self.target.chains) * [dc]
    else:
      self.dc = None

    #set all the other parameters
    self.potential_container=potential_container
    self.smooth_std=smooth_std
    self.norm=norm
    self.ca_positions=None
    self.spherical_smoother=None
    self.count_radius = count_radius
    
    if target.handle != environment.handle:
      raise RuntimeError('Target and Environment must be views of the same handle!')

    #if psipred data is provided, take secondary information from there. Particularly in models with
    #low quality, the secondary structure assignment of DSSP srews up. In case of membrane proteins,
    #the secondary structure of the membrane spanning region is either set to helical or extended, 
    #depending on the majority of residues.
    mol.alg.AssignSecStruct(environment)
    mol.alg.Accessibility(environment, algorithm=mol.alg.DSSP)

    dssp_ss = list()

    for r in target.residues:
      if r.GetSecStructure().IsHelical():
        dssp_ss+='H'
      elif r.GetSecStructure().IsExtended():
        dssp_ss+='E'
      else:
        dssp_ss+='C' 

    if psipred != None:
      self.sec_structure = ''
      for c,p in zip(self.target.chains, self.psipred):
        if p != None:
          self.sec_structure+= p.GetPSIPREDSS(c)
        else:
          self.sec_structure+=''.join(str(r.GetSecStructure()) for r in c.residues)
    else:
      self.sec_structure = dssp_ss
 
    #map the secondary structure onto the target
    for ss, r in zip(self.sec_structure, self.target.residues):
      if ss=='H':
        r.SetIntProp('SS',0)
      elif ss=='E':
        r.SetIntProp('SS',1)
      else:
        r.SetIntProp('SS',2)

    #get ca positions, that are used for spherical smoothing later on
    if smooth_std!=None:
      pos=list()
      for r in self.target.residues:
        ca = r.FindAtom('CA')
        if ca.IsValid():
          pos.append(ca.GetPos())
        else:
          pos.append(r.GetCenterOfMass())
      self.ca_positions=pos
      self.spherical_smoother=SphericalSmoother(self.ca_positions,smooth_std)

    #make the selections...
    self.helix_selection = self.target.Select('grSS=0')
    self.extended_selection = self.target.Select('grSS=1')
    self.coil_selection = self.target.Select('grSS=2')

    sequence=list()
    for r in self.target.residues:
      sequence.append(r.one_letter_code)

    self.data['sequence'] = sequence
    self.data['sec_structure'] = self.sec_structure
    self.data['dssp_ss'] = dssp_ss

  def GetLocalData(self, features):

    return_data = dict()
    for f in features:
      if f not in self.data:
        self.Calculate(f)
      return_data[f] = self.data[f]
    return return_data

  def GetAVGData(self, features):

    return_data = dict()
    for f in features:
      if f not in self.data:
        self.Calculate(f)
      return_data[f] = self.data["avg_" + f]
    return return_data

  def Calculate(self, feature):
    if feature == 'torsion':
      self.DoTorsion()
    elif feature == 'interaction':
      self.DoInteraction()
    elif feature == 'cbeta':
      self.DoCBeta()
    elif feature == 'reduced':
      self.DoReduced()
    elif feature == 'packing':
      self.DoPacking()
    elif feature == 'cb_packing':
      self.DoCBPacking()
    elif feature == 'ss_agreement':
      self.DoSSAgreement()
    elif feature == 'acc_agreement':
      self.DoACCAgreement()
    elif feature == 'dist_const':
      self.DoConstraints()
    elif feature == 'fraction_loops':
      self.DoFractionLoops()
    elif feature == 'exposed':
      self.DoExposed()
    elif feature == 'counts':
      self.DoCounts()
    elif feature == 'clash':
      self.DoClash()
    else:
      raise ValueError('Requested feature \"%s\" is not supported!' %(f))

  def GetTarget(self):
    return self.target

  def GetEnvironment(self):
    return self.environment

  def GetDSSPSS(self):
    return self.data['dssp_ss']

  def GetAssignedSS(self):
    return self.data['sec_structure']

  def DoTorsion(self):
    torsion = self.potential_container['torsion'].GetEnergies(self.target)
    self.data['avg_torsion'] = self.GetAverage(torsion)
    if self.smooth_std!=None:
      self.data['torsion'] = self.spherical_smoother.Smooth(torsion)
    else:
      self.data['torsion'] = torsion

  def DoInteraction(self, ss_dependent = True):

    interaction = None

    if ss_dependent:
      helix_interaction = self.potential_container['interaction_helix'].GetEnergies(self.helix_selection, self.environment, normalize=self.norm)
      extended_interaction = self.potential_container['interaction_extended'].GetEnergies(self.extended_selection, self.environment, normalize=self.norm)
      coil_interaction = self.potential_container['interaction_coil'].GetEnergies(self.coil_selection, self.environment, normalize=self.norm)

      interaction = list()

      for r in self.target.residues:
        if r.GetIntProp('SS')==0:
          interaction.append(helix_interaction.pop(0))
        elif r.GetIntProp('SS')==1:
          interaction.append(extended_interaction.pop(0))
        else:
          interaction.append(coil_interaction.pop(0))
    else:
      interaction = self.potential_container['interaction'].GetEnergies(self.target, self.environment, normalize = self.norm)

    self.data['avg_interaction'] = self.GetAverage(interaction)

    if self.smooth_std!=None:
      self.data['interaction'] = self.spherical_smoother.Smooth(interaction)
    else:
      self.data['interaction'] = interaction

  def DoCBeta(self, ss_dependent = True):

    cbeta = None

    if ss_dependent:
      helix_cbeta = self.potential_container['cbeta_helix'].GetEnergies(self.helix_selection, self.environment, normalize=self.norm)
      extended_cbeta = self.potential_container['cbeta_extended'].GetEnergies(self.extended_selection, self.environment, normalize=self.norm)
      coil_cbeta = self.potential_container['cbeta_coil'].GetEnergies(self.coil_selection, self.environment, normalize=self.norm)

      cbeta = list()

      for r in self.target.residues:
        if r.GetIntProp('SS')==0:
          cbeta.append(helix_cbeta.pop(0))
        elif r.GetIntProp('SS')==1:
          cbeta.append(extended_cbeta.pop(0))
        else:
          cbeta.append(coil_cbeta.pop(0))

    else:
      cbeta = self.potential_container['cbeta'].GetEnergies(self.target, self.environment, normalize = self.norm)

    self.data['avg_cbeta'] = self.GetAverage(cbeta)

    if self.smooth_std!=None:
      self.data['cbeta'] = self.spherical_smoother.Smooth(cbeta)
    else:
      self.data['cbeta'] = cbeta

  def DoReduced(self, ss_dependent = True):

    reduced = None

    if ss_dependent:
      helix_reduced = self.potential_container['reduced_helix'].GetEnergies(self.helix_selection, self.environment, normalize=self.norm)
      extended_reduced = self.potential_container['reduced_extended'].GetEnergies(self.extended_selection, self.environment, normalize=self.norm)
      coil_reduced = self.potential_container['reduced_coil'].GetEnergies(self.coil_selection, self.environment, normalize=self.norm)

      reduced = list()

      for r in self.target.residues:
        if r.GetIntProp('SS')==0:
          reduced.append(helix_reduced.pop(0))
        elif r.GetIntProp('SS')==1:
          reduced.append(extended_reduced.pop(0))
        else:
          reduced.append(coil_reduced.pop(0))
    else:
      reduced = self.potential_container['reduced'].GetEnergies(self.target, self.environment, normalize = self.norm)
      

    self.data['avg_reduced'] = self.GetAverage(reduced)

    if self.smooth_std!=None:
      self.data['reduced'] = self.spherical_smoother.Smooth(reduced)
    else:
      self.data['reduced'] = reduced

  def DoPacking(self, ss_dependent = True):

    packing=None

    if ss_dependent:
      helix_packing = self.potential_container['packing_helix'].GetEnergies(self.helix_selection, self.environment,normalize=self.norm)
      extended_packing = self.potential_container['packing_extended'].GetEnergies(self.extended_selection, self.environment,normalize=self.norm)
      coil_packing = self.potential_container['packing_coil'].GetEnergies(self.coil_selection, self.environment,normalize=self.norm)

      packing = list()

      for r in self.target.residues:
        if r.GetIntProp('SS')==0:
          packing.append(helix_packing.pop(0))
        elif r.GetIntProp('SS')==1:
          packing.append(extended_packing.pop(0))
        else:
          packing.append(coil_packing.pop(0))
    else:
      packing = self.potential_container['packing'].GetEnergies(self.target, self.environment, normalize = self.norm)


    self.data['avg_packing'] = self.GetAverage(packing)

    if self.smooth_std!=None:
      self.data['packing'] = self.spherical_smoother.Smooth(packing)
    else:
      self.data['packing'] = packing

  def DoCBPacking(self, ss_dependent = True):

    cb_packing=None

    if ss_dependent:
      helix_cb_packing = self.potential_container['cb_packing_helix'].GetEnergies(self.helix_selection, self.environment)
      extended_cb_packing = self.potential_container['cb_packing_extended'].GetEnergies(self.extended_selection, self.environment)
      coil_cb_packing = self.potential_container['cb_packing_coil'].GetEnergies(self.coil_selection, self.environment)

      cb_packing = list()

      for r in self.target.residues:
        if r.GetIntProp('SS')==0:
          cb_packing.append(helix_cb_packing.pop(0))
        elif r.GetIntProp('SS')==1:
          cb_packing.append(extended_cb_packing.pop(0))
        else:
          cb_packing.append(coil_cb_packing.pop(0))
    else:
      cb_packing = self.potential_container['cb_packing'].GetEnergies(self.target, self.environment)

    self.data['avg_cb_packing'] = self.GetAverage(cb_packing)

    if self.smooth_std!=None:
      self.data['cb_packing'] = self.spherical_smoother.Smooth(cb_packing)
    else:
      self.data['cb_packing'] = cb_packing

  def DoExposed(self):
    exposed=list()
    for r in self.target.residues:
      try:
        exposed.append(r.GetFloatProp('asaRel'))
      except:
        exposed.append(float('NaN'))
    self.data['avg_exposed'] = self.GetAverage(exposed)
    if self.smooth_std!=None:
      self.data['exposed']=self.spherical_smoother.Smooth(exposed)
    else:
      self.data['exposed']=exposed

  def DoFractionLoops(self):
    fraction_loops=list()

    for r in self.target.residues:
      if r.GetSecStructure().IsCoil():
        fraction_loops.append(1.0)
      else:
        fraction_loops.append(0.0)

    self.data['avg_fraction_loops'] = self.GetAverage(fraction_loops)
    if self.smooth_std!=None:
      self.data['fraction_loops']=self.spherical_smoother.Smooth(fraction_loops)
    else:
      self.data['fraction_loops']=fraction_loops

  def DoCounts(self):
    counts=list()

    ca_selection = self.environment.Select("aname=CA")
    for r in self.target.residues:
      ca = r.FindAtom("CA")
      if ca.IsValid():
        close_stuff = ca_selection.FindWithin(ca.GetPos(), self.count_radius)
        counts.append(float(len(close_stuff) - 1))
      else:
        counts.append(float("NaN"))
    self.data['avg_counts'] = self.GetAverage(counts)
    self.data['counts'] = counts

  def DoSSAgreement(self):
    if self.psipred!=None:
      ss_agreement = []
      for c, p in zip(self.target.chains, self.psipred):
        ss_agreement += p.GetSSAgreementFromChain(c)
      self.data['avg_ss_agreement'] = self.GetAverage(ss_agreement)
      if self.smooth_std!=None:
        self.data['ss_agreement']=self.spherical_smoother.Smooth(ss_agreement)
      else:
        self.data['ss_agreement']=ss_agreement
    else:
      self.data["avg_ss_agreement"] = float("NaN")
      self.data["ss_agreement"] = [float("NaN")] * len(self.target.residues)

  def DoACCAgreement(self):
    if self.accpro!=None:
      acc_agreement = []
      for c,a in zip(self.target.chains, self.accpro):
        acc_agreement += a.GetACCAgreementFromChain(c)
      self.data['avg_acc_agreement'] = self.GetAverage(acc_agreement)
      if self.smooth_std!=None:
        self.data['acc_agreement']=self.spherical_smoother.Smooth(acc_agreement)
      else:
        self.data['acc_agreement']=acc_agreement
    else:
      self.data["avg_acc_agreement"] = float("NaN")
      self.data["acc_agreement"] = [float("NaN")] * len(self.target.residues)

  def DoClash(self):
    clash_scores = GetClashScores(self.target, self.environment)
    self.data["clash"] = clash_scores
    self.data["avg_clash"] = self.GetAverage(clash_scores)
    

  def DoConstraints(self):
    if self.dc!=None:

      dist_const_data = dict()
      dist_const_data["disco"] = list()
      dist_const_data["counts"] = list()
      dist_const_data["avg_num_clusters"] = list()
      dist_const_data["avg_max_seqsim"] = list()
      dist_const_data["avg_max_seqid"] = list()
      dist_const_data["avg_variance"] = list()
      dist_const_data["num_constraints"] = list()
      dist_const_data["fraction_observed"] = list()

      for c,d in zip(self.target.chains, self.dc):

        if d!= None:

          entity_view = self.target.CreateEmptyView()
          entity_view.AddChain(c, mol.INCLUDE_ALL)
          chain_data = d.GetScoresWithData(entity_view)
          dist_const_data["disco"] += chain_data["scores"]
          dist_const_data["counts"] += chain_data["counts"]
          dist_const_data["avg_num_clusters"] += chain_data["avg_num_clusters"]
          dist_const_data["avg_max_seqsim"] += chain_data["avg_max_seqsim"]
          dist_const_data["avg_max_seqid"] += chain_data["avg_max_seqid"]
          dist_const_data["avg_variance"] += chain_data["avg_variance"]
          dist_const_data["num_constraints"] += chain_data["num_constraints"]
 
          fraction_observed = list()
          for a,b in zip(chain_data["counts"], chain_data["num_constraints"]):
            if b == 0:
              fraction_observed.append(0.0)
            else:
              fraction_observed.append(float(a)/b)
          dist_const_data["fraction_observed"] += fraction_observed

        else:

          dist_const_data["disco"] += list([float("NaN")] * len(c.residues))
          dist_const_data["counts"] += list([float("NaN")] * len(c.residues))
          dist_const_data["avg_num_clusters"] += list([float("NaN")] * len(c.residues))
          dist_const_data["avg_max_seqsim"] += list([float("NaN")] * len(c.residues))
          dist_const_data["avg_max_seqid"] += list([float("NaN")] * len(c.residues))
          dist_const_data["avg_variance"] += list([float("NaN")] * len(c.residues))
          dist_const_data["num_constraints"] += list([float("NaN")] * len(c.residues))
          dist_const_data["fraction_observed"] += list([float("NaN")] * len(c.residues))

      self.data["dist_const"] = dist_const_data
      self.data["avg_dist_const"] = self.GetAverage(dist_const_data["disco"])   

    else:   
      dist_const_data = dict()
      nan_list = [float("NaN")] * len(self.target.residues)
      dist_const_data["disco"] = nan_list
      dist_const_data["counts"] = nan_list
      dist_const_data["avg_num_clusters"] = nan_list
      dist_const_data["avg_max_seqsim"] = nan_list
      dist_const_data["avg_max_seqid"] = nan_list
      dist_const_data["avg_variance"] = nan_list
      dist_const_data["num_constraints"] = nan_list
      self.data["dist_const"] = dist_const_data
      self.data["avg_dist_const"] = float("NaN")
    
  def GetAverage(self, value_list):
    temp = list()
    for item in value_list:
      if item==item:
        temp.append(float(item))
    if len(temp)>0:
      return sum(temp)/len(temp)
    else:
      return float('NaN')
