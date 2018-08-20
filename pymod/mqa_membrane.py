# QMEAN:
# ----------
# 
# Copyright (C) 2018 University of Basel and the Authors.
# Authors: Gabriel Studer, Pascal Benkert, Marco Biasini
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

from qmean import *
from qmean import conf
from ost import mol
import pickle
from predicted_sequence_features import AlignChainToSEQRES


class MembraneScores:

  def __init__(self, target, environment, potential_container_soluble, potential_container_membrane, smooth_std=None, psipred=None, 
               accpro=None, norm=True, mem_param=None, membrane_query = None, interface_query=None):
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

    if accpro != None:
      if isinstance(accpro,list):
        if len(accpro) != len(target.chains):
          raise ValueError("List of provided ACCPROHandlers must have same length as chains in target!")
        self.accpro = accpro
      else:
        self.accpro = len(self.target.chains) * [accpro]

    #set all the other parameters
    self.potential_container_soluble=potential_container_soluble
    self.potential_container_membrane=potential_container_membrane
    self.smooth_std=smooth_std
    self.norm=norm
    self.mem_param = mem_param
    self.ca_positions=None
    self.spherical_smoother=None

    for r in self.target.residues:
      r.SetIntProp('membrane_state',0)
      r.SetIntProp('in_membrane',0)

    if membrane_query != None and interface_query != None:
      membrane_selection = target.Select(membrane_query)
      interface_selection = target.Select(interface_query)
      intersection = mol.Intersection(membrane_selection, interface_selection)
      if len(intersection.residues) > 0:
        raise RuntimeError("The views resulting from your membrane and interface query must not intersect!")
      for r in membrane_selection.residues:
        r.SetIntProp('membrane_state',1)
      for r in interface_selection.residues:
        r.SetIntProp('membrane_state',2)
      
      ca_selection = self.target.Select('aname=CA')

      for r in membrane_selection.residues:
        ca = r.FindAtom('CA')
        if ca.IsValid():
          in_range = ca_selection.FindWithin(ca.GetPos(),5.0)
          for a in in_range:
            a.GetResidue().SetIntProp('in_membrane',1)

    else:
      if self.mem_param == None:
        self.mem_param = mol.alg.FindMembrane(self.target)

      membrane_axis = geom.Normalize(self.mem_param.membrane_axis)
      membrane_center = self.mem_param.pos
      membrane_width = self.mem_param.width

      top = membrane_center * membrane_axis + membrane_axis * membrane_width/2
      bottom = membrane_center * membrane_axis - membrane_axis * membrane_width/2

      plane_one = geom.Plane(top, membrane_axis)
      plane_two = geom.Plane(bottom, membrane_axis)

      for r in self.target.residues:
        ca = r.FindAtom('CA')
        if not ca.IsValid():
          continue
        dist_one = abs(geom.Distance(plane_one,ca.GetPos()))
        dist_two = abs(geom.Distance(plane_two,ca.GetPos()))
        if dist_one < membrane_width and dist_two < membrane_width:
          r.SetIntProp('in_membrane', 1)
        if dist_one < 5 or dist_two < 5:
          r.SetIntProp('membrane_state',2)
          continue
        if dist_one < membrane_width and dist_two < membrane_width:
          r.SetIntProp('membrane_state',1)

    self.interface_target = self.target.Select('grmembrane_state=2')
    self.membrane_target = self.target.Select('grmembrane_state=1')
    self.soluble_target = self.target.Select('grmembrane_state=0')

    self.membrane_states = list()
    self.membrane_associations = list()

    for r in self.target.residues:
      self.membrane_states.append(r.GetIntProp('membrane_state'))
      self.membrane_associations.append(r.GetIntProp('in_membrane'))

    if target.handle != environment.handle:
      raise RuntimeError('Target and Environment must be views of the same handle!')

    mol.alg.AssignSecStruct(target)
    mol.alg.Accessibility(target, algorithm=mol.alg.DSSP)

    #if there is a membrane, check whether the structure has a helical or extended
    #transmembranepart

    num_helical_membrane_residues=0
    num_extended_membrane_residues=0

    for r in self.membrane_target.residues:
      if r.GetSecStructure().IsHelical():
        num_helical_membrane_residues+=1
      elif r.GetSecStructure().IsExtended():
        num_extended_membrane_residues+=1

    if num_helical_membrane_residues < num_extended_membrane_residues:
      raise RuntimeError("QMEANBrane is currently only trained on alpha helical membrane proteins!")

    #if psipred data is provided, take secondary information from there. Particularly in models with
    #low quality, the secondary structure assignment of DSSP screws up. In case of membrane proteins,
    #the secondary structure of the membrane spanning region is either set to helical or extended, 
    #depending on the majority of residues.

    self.dssp_ss = list()

    for r in target.residues:
      if r.GetSecStructure().IsHelical():
        self.dssp_ss.append('H')
      elif r.GetSecStructure().IsExtended():
        self.dssp_ss.append('E')
      else:
        self.dssp_ss.append('C')

    temp_ss = ""
    if psipred != None:
      for c,p in zip(self.target.chains,self.psipred):
        temp_ss += p.GetPSIPREDSS(c)
    else:
      temp_ss = self.dssp_ss

    self.sec_structure=list()
    for ss, r in zip(temp_ss, self.target.residues):
      if self.membrane_target.ViewForHandle(r.handle).IsValid():
        self.sec_structure.append('H')
      else:
        self.sec_structure.append(ss)

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
    self.soluble_helix_target = self.soluble_target.Select('grSS=0')
    self.soluble_extended_target = self.soluble_target.Select('grSS=1')
    self.soluble_coil_target = self.soluble_target.Select('grSS=2')

    sequence=list()
    for r in self.target.residues:
      sequence.append(r.one_letter_code)

    self.data['sequence'] = sequence
    self.data['sec_structure'] = self.sec_structure
    self.data['dssp_ss'] = self.dssp_ss


  def GetLocalData(self,features):
    return_data = dict()
    for f in features:
      if f == 'torsion':
        try:
          return_data['torsion']=self.data['torsion']
        except:
          raise ValueError("Did not calculate the torsion feature!")
      elif f == 'interaction':
        try:
          return_data['interaction']=self.data['interaction']
        except:
          raise ValueError("Did not calculate the interaction feature!")
      elif f == 'cbeta':
        try:
          return_data['cbeta']=self.data['cbeta']
        except:
          raise ValueError("Did not calculate the cbeta feature!")
      elif f == 'packing':
        try:
          return_data['packing']=self.data['packing']
        except:
          raise ValueError("Did not calculate the packing feature!")
      elif f == 'ss_agreement':
        try:
          return_data['ss_agreement'] = self.data['ss_agreement']
        except:
          raise ValueError("Did not calculate the ss_agreement feature!")
      elif f == 'acc_agreement':
        try:
          return_data['acc_agreement'] = self.data['acc_agreement']
        except:
          raise ValueError("Did not calculate the acc_agreement feature!")   
      elif f == 'fraction_loops':
        try:
          return_data['fraction_loops'] = self.data['fraction_loops']
        except:
          raise ValueError("Did not calculate the fraction_loops feature!")
      elif f == 'exposed':
        try:
          return_data['exposed'] = self.data['exposed']
        except:
          raise ValueError("Did not calculate the exposed feature!")
      else:
        raise ValueError('Requested feature is not supported!')
        
    return return_data

  def GetAVGData(self,features):
    return_data = dict()
    for f in features:
      if f == 'torsion':
        try:
          return_data['torsion']=self.data['avg_torsion']
        except:
          raise ValueError("Did not calculate the torsion feature!")
      elif f == 'interaction':
        try:
          return_data['interaction']=self.data['avg_interaction']
        except:
          raise ValueError("Did not calculate the interaction feature!")
      elif f == 'cbeta':
        try:
          return_data['cbeta']=self.data['avg_cbeta']
        except:
          raise ValueError("Did not calculate the cbeta feature!")
      elif f == 'packing':
        try:
          return_data['packing']=self.data['avg_packing']
        except:
          raise ValueError("Did not calculate the packing feature!")
      elif f == 'ss_agreement':
        try:
          return_data['ss_agreement'] = self.data['avg_ss_agreement']
        except:
          raise ValueError("Did not calculate the ss_agreement feature!")
      elif f == 'acc_agreement':
        try:
          return_data['acc_agreement'] = self.data['avg_acc_agreement']
        except:
          raise ValueError("Did not calculate the acc_agreement feature!")
      elif f == 'fraction_loops':
        try:
          return_data['fraction_loops'] = self.data['avg_fraction_loops']
        except:
          raise ValueError("Did not calculate the fraction_loops feature!")
      elif f == 'exposed':
        try:
          return_data['exposed'] = self.data['avg_exposed']
        except:
          raise ValueError("Did not calculate the exposed feature!")
      else:
        raise ValueError('Requested feature \"%s\" is not supported!' %(f))
        
    return return_data

  def GetTarget(self):
    return self.target

  def GetEnvironment(self):
    return self.environment

  def GetSolubleTarget(self):
    return self.soluble_target

  def GetMembraneTarget(self):
    return self.membrane_target

  def GetInterfaceTarget(self):
    return self.interface_target

  def GetDSSPSS(self):
    return self.dssp_ss

  def GetAssignedSS(self):
    return self.sec_structure

  def GetMembraneState(self):
    return self.membrane_states

  def GetMembraneAssociation(self):
    return self.membrane_associations

  def GetMemParam(self):
    return self.mem_param

  def CalculateScores(self, features):

    for f in features:
      if f == 'torsion':
        self.GetTorsion()
      elif f == 'interaction':
        self.GetInteraction()
      elif f == 'cbeta':
        self.GetCBeta()
      elif f == 'packing':
        self.GetPacking()
      elif f == 'ss_agreement':
        if self.psipred!=None:
          self.GetSSAgreement()
        else:
          raise ValueError("you must provide psipred data for the secondary structure agreement term!")
      elif f == 'acc_agreement':
        if self.accpro!=None:
          self.GetACCAgreement()
        else:
          raise ValueError("you must provide accpro data for acc_agreement term!")  
      elif f == 'fraction_loops':
        self.GetFractionLoops()
      elif f == 'exposed':
        self.GetExposed()
      else:
        raise ValueError('Requested feature is not supported in local scorer!')

  def GetTorsion(self):
    #we use the soluble torsion potential in the membrane!
    torsion = self.potential_container_soluble['torsion'].GetEnergies(self.target)
    self.data['avg_torsion'] = self.GetAverage(torsion)
    if self.smooth_std!=None:
      self.data['torsion'] = self.spherical_smoother.Smooth(torsion)
    else:
      self.data['torsion'] = torsion
 
  def GetInteraction(self):

    interaction=list()

    soluble_helix_interaction = self.potential_container_soluble['interaction_helix'].GetEnergies(self.soluble_helix_target, self.environment, normalize=self.norm)
    soluble_extended_interaction = self.potential_container_soluble['interaction_extended'].GetEnergies(self.soluble_extended_target, self.environment, normalize=self.norm)
    soluble_coil_interaction = self.potential_container_soluble['interaction_coil'].GetEnergies(self.soluble_coil_target, self.environment, normalize=self.norm)

    membrane_interaction = self.potential_container_membrane['membrane_interaction_helix'].GetEnergies(self.membrane_target,self.environment,normalize=self.norm)
    interface_interaction = self.potential_container_membrane['interface_interaction_helix'].GetEnergies(self.interface_target,self.environment,normalize=self.norm)

    for r in self.target.residues:
      if self.membrane_target.ViewForHandle(r.handle).IsValid():
        interaction.append(membrane_interaction.pop(0))
      elif self.interface_target.ViewForHandle(r.handle).IsValid():
        interaction.append(interface_interaction.pop(0))
      else:
        if r.GetIntProp('SS')==0:
          interaction.append(soluble_helix_interaction.pop(0))
        elif r.GetIntProp('SS')==1:
          interaction.append(soluble_extended_interaction.pop(0))
        else:
          interaction.append(soluble_coil_interaction.pop(0))

    if self.smooth_std!=None:
      self.data['interaction'] = self.spherical_smoother.Smooth(interaction)
    else:
      self.data['interaction'] = interaction

    self.data['avg_interaction'] = self.GetAverage(interaction)


  def GetCBeta(self):
    cbeta=list()

    soluble_helix_cbeta = self.potential_container_soluble['cbeta_helix'].GetEnergies(self.soluble_helix_target, self.environment, normalize=self.norm)
    soluble_extended_cbeta = self.potential_container_soluble['cbeta_extended'].GetEnergies(self.soluble_extended_target, self.environment, normalize=self.norm)
    soluble_coil_cbeta = self.potential_container_soluble['cbeta_coil'].GetEnergies(self.soluble_coil_target, self.environment, normalize=self.norm)

    membrane_cbeta = self.potential_container_membrane['membrane_cbeta_helix'].GetEnergies(self.membrane_target,self.environment,normalize=self.norm)
    interface_cbeta = self.potential_container_membrane['interface_cbeta_helix'].GetEnergies(self.interface_target,self.environment,normalize=self.norm)

    for r in self.target.residues:
      if self.membrane_target.ViewForHandle(r.handle).IsValid():
        cbeta.append(membrane_cbeta.pop(0))
      elif self.interface_target.ViewForHandle(r.handle).IsValid():
        cbeta.append(interface_cbeta.pop(0))
      else:
        if r.GetIntProp('SS')==0:
          cbeta.append(soluble_helix_cbeta.pop(0))
        elif r.GetIntProp('SS')==1:
          cbeta.append(soluble_extended_cbeta.pop(0))
        else:
          cbeta.append(soluble_coil_cbeta.pop(0))

    if self.smooth_std!=None:
      self.data['cbeta'] = self.spherical_smoother.Smooth(cbeta)
    else:
      self.data['cbeta'] = cbeta

    self.data['avg_cbeta'] = self.GetAverage(cbeta)


  def GetPacking(self):
    packing=list()

    soluble_helix_packing = self.potential_container_soluble['packing_helix'].GetEnergies(self.soluble_helix_target, self.environment,normalize=self.norm)
    soluble_extended_packing = self.potential_container_soluble['packing_extended'].GetEnergies(self.soluble_extended_target, self.environment,normalize=self.norm)
    soluble_coil_packing = self.potential_container_soluble['packing_coil'].GetEnergies(self.soluble_coil_target, self.environment,normalize=self.norm)

    membrane_packing = self.potential_container_membrane['membrane_packing_helix'].GetEnergies(self.membrane_target,self.environment,normalize=self.norm)
    interface_packing = self.potential_container_membrane['interface_packing_helix'].GetEnergies(self.interface_target,self.environment,normalize=self.norm)

    for r in self.target.residues:
      if self.membrane_target.ViewForHandle(r.handle).IsValid():
        packing.append(membrane_packing.pop(0))
      elif self.interface_target.ViewForHandle(r.handle).IsValid():
        packing.append(interface_packing.pop(0))
      else:
        if r.GetIntProp('SS')==0:
          packing.append(soluble_helix_packing.pop(0))
        elif r.GetIntProp('SS')==1:
          packing.append(soluble_extended_packing.pop(0))
        else:
          packing.append(soluble_coil_packing.pop(0))

    if self.smooth_std!=None:
      self.data['packing'] = self.spherical_smoother.Smooth(packing)
    else:
      self.data['packing'] = packing

    self.data['avg_packing'] = self.GetAverage(packing)

  
  def GetExposed(self):

    exposed=list()
    for r in self.target.residues:
      try:
        exposed.append(r.GetFloatProp('asaRel'))
      except:
        exposed.append(float('NaN'))
    if self.smooth_std!=None:
      self.data['exposed']=self.spherical_smoother.Smooth(exposed)
    else:
      self.data['exposed']=exposed
    self.data['avg_exposed'] = self.GetAverage(exposed)

  def GetFractionLoops(self):
    fraction_loops=list()
    for r in self.target.residues:
      if r.GetSecStructure().IsCoil():
        fraction_loops.append(1.0)
      else:
        fraction_loops.append(0.0)
    if self.smooth_std!=None:
      self.data['fraction_loops']=self.spherical_smoother.Smooth(fraction_loops)
    else:
      self.data['fraction_loops']=fraction_loops

    self.data['avg_fraction_loops'] = self.GetAverage(fraction_loops)

  def GetSSAgreement(self):
    if self.psipred!=None:
      ss_agreement = []
      for c,p in zip(self.target.chains,self.psipred):
        ss_agreement += p.GetSSAgreementFromChain(c)
      self.data['avg_ss_agreement'] = self.GetAverage(ss_agreement)
      if self.smooth_std!=None:
        self.data['ss_agreement']=self.spherical_smoother.Smooth(ss_agreement)
      else:
        self.data['ss_agreement']=ss_agreement
    else:
      raise ValueError("Cannot calculate ss_agreement term without psipred information!")

  def GetACCAgreement(self):
    if self.accpro!=None:
      acc_agreement = []
      for c,a in zip(self.target.chains,self.accpro):
        acc_agreement += a.GetACCAgreementFromChain(c)
      self.data['avg_acc_agreement'] = self.GetAverage(acc_agreement)
      if self.smooth_std!=None:
        self.data['acc_agreement']=self.spherical_smoother.Smooth(acc_agreement)
      else:
        self.data['acc_agreement']=acc_agreement
    else:
      raise ValueError("Cannot calculate acc_agreement term without accpro information!")

  def GetAverage(self, value_list):
    temp = list()
    for item in value_list:
      if item==item:
        temp.append(float(item))
    if len(temp)>0:
      return sum(temp)/len(temp)
    else:
      return float('NaN')
