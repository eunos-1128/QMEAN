from qmean import *
from ost.bindings import dssp
from ost import mol


class Scores:

  def __init__(self, target, environment, potential_container, smooth_std=None, psipred=None, accpro=None, assign_dssp=True, norm=True):
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
    self.potential_container=potential_container
    self.smooth_std=smooth_std
    self.norm=norm
    self.ca_positions=None
    self.spherical_smoother=None
    
    if target.handle != environment.handle:
      raise RuntimeError('Target and Environment must be views of the same handle!')

    if assign_dssp:
      dssp.AssignDSSP(environment, extract_burial_status=True)

    #if psipred data is provided, take secondary information from there. Particularly in models with
    #low quality, the secondary structure assignment of DSSP srews up. In case of membrane proteins,
    #the secondary structure of the membrane spanning region is either set to helical or extended, 
    #depending on the majority of residues.

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
        self.sec_structure+= p.GetPSIPREDSS(c)
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
      elif f == 'reduced':
        try:
          return_data['reduced']=self.data['reduced']
        except:
          raise ValueError("Did not calculate the reduced feature!")
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
      elif f == 'reduced':
        try:
          return_data['reduced']=self.data['avg_reduced']
        except:
          raise ValueError("Did not calculate the reduced feature!")
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

  def GetDSSPSS(self):
    return self.data['dssp_ss']

  def GetAssignedSS(self):
    return self.data['sec_structure']

  def CalculateScores(self, features):
    for f in features:
      if f == 'torsion':
        self.GetTorsion()
      elif f == 'interaction':
        self.GetInteraction()
      elif f == 'cbeta':
        self.GetCBeta()
      elif f == 'reduced':
        self.GetReduced()
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
          raise ValueError("you must provide accpro data for ACC Agree term!")
      elif f == 'fraction_loops':
        self.GetFractionLoops()
      elif f == 'exposed':
        self.GetExposed()
      else:
        raise ValueError('Requested feature \"%s\" is not supported!' %(f))

  def GetTorsion(self):
    torsion = self.potential_container['torsion'].GetEnergies(self.target, normalize=self.norm)
    self.data['avg_torsion'] = self.GetAverage(torsion)
    if self.smooth_std!=None:
      self.data['torsion'] = self.spherical_smoother.Smooth(torsion)
    else:
      self.data['torsion'] = torsion

  def GetInteraction(self):
    interaction=list()

    helix_interaction = self.potential_container['interaction_helix'].GetEnergies(self.helix_selection, self.environment, normalize=self.norm)
    extended_interaction = self.potential_container['interaction_extended'].GetEnergies(self.extended_selection, self.environment, normalize=self.norm)
    coil_interaction = self.potential_container['interaction_coil'].GetEnergies(self.coil_selection, self.environment, normalize=self.norm)

    for r in self.target.residues:
      if r.GetIntProp('SS')==0:
        interaction.append(helix_interaction.pop(0))
      elif r.GetIntProp('SS')==1:
        interaction.append(extended_interaction.pop(0))
      else:
        interaction.append(coil_interaction.pop(0))

    self.data['avg_interaction'] = self.GetAverage(interaction)

    if self.smooth_std!=None:
      self.data['interaction'] = self.spherical_smoother.Smooth(interaction)
    else:
      self.data['interaction'] = interaction

  def GetCBeta(self):
    cbeta=list()

    helix_cbeta = self.potential_container['cbeta_helix'].GetEnergies(self.helix_selection, self.environment, normalize=self.norm)
    extended_cbeta = self.potential_container['cbeta_extended'].GetEnergies(self.extended_selection, self.environment, normalize=self.norm)
    coil_cbeta = self.potential_container['cbeta_coil'].GetEnergies(self.coil_selection, self.environment, normalize=self.norm)

    for r in self.target.residues:
      if r.GetIntProp('SS')==0:
        cbeta.append(helix_cbeta.pop(0))
      elif r.GetIntProp('SS')==1:
        cbeta.append(extended_cbeta.pop(0))
      else:
        cbeta.append(coil_cbeta.pop(0))

    self.data['avg_cbeta'] = self.GetAverage(cbeta)

    if self.smooth_std!=None:
      self.data['cbeta'] = self.spherical_smoother.Smooth(cbeta)
    else:
      self.data['cbeta'] = cbeta

  def GetReduced(self):
    reduced=list()

    helix_reduced = self.potential_container['reduced_helix'].GetEnergies(self.helix_selection, self.environment, normalize=self.norm)
    extended_reduced = self.potential_container['reduced_extended'].GetEnergies(self.extended_selection, self.environment, normalize=self.norm)
    coil_reduced = self.potential_container['reduced_coil'].GetEnergies(self.coil_selection, self.environment, normalize=self.norm)

    for r in self.target.residues:
      if r.GetIntProp('SS')==0:
        reduced.append(helix_reduced.pop(0))
      elif r.GetIntProp('SS')==1:
        reduced.append(extended_reduced.pop(0))
      else:
        reduced.append(coil_reduced.pop(0))

    self.data['avg_reduced'] = self.GetAverage(reduced)

    if self.smooth_std!=None:
      self.data['reduced'] = self.spherical_smoother.Smooth(reduced)
    else:
      self.data['reduced'] = reduced

  def GetPacking(self):
    packing=list()

    helix_packing = self.potential_container['packing_helix'].GetEnergies(self.helix_selection, self.environment,normalize=self.norm)
    extended_packing = self.potential_container['packing_extended'].GetEnergies(self.extended_selection, self.environment,normalize=self.norm)
    coil_packing = self.potential_container['packing_coil'].GetEnergies(self.coil_selection, self.environment,normalize=self.norm)

    for r in self.target.residues:
      if r.GetIntProp('SS')==0:
        packing.append(helix_packing.pop(0))
      elif r.GetIntProp('SS')==1:
        packing.append(extended_packing.pop(0))
      else:
        packing.append(coil_packing.pop(0))

    self.data['avg_packing'] = self.GetAverage(packing)

    if self.smooth_std!=None:
      self.data['packing'] = self.spherical_smoother.Smooth(packing)
    else:
      self.data['packing'] = packing

  def GetExposed(self):
    exposed=list()

    for r in self.target.residues:
      try:
        exposed.append(r.GetFloatProp('relative_solvent_accessibility'))
      except:
        exposed.append(float('NaN'))

    self.data['avg_exposed'] = self.GetAverage(exposed)

    if self.smooth_std!=None:
      self.data['exposed']=self.spherical_smoother.Smooth(exposed)
    else:
      self.data['exposed']=exposed

  def GetFractionLoops(self):
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

  def GetSSAgreement(self):
    if self.psipred!=None:
      ss_agreement = []
      for c, p in zip(self.target.chains, self.psipred):
        ss_agreement += p.GetSSAgreementFromChain(c, dssp_assigned=True)
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
      for c,a in in zip(self.target.chains, self.accpro):
        acc_agreement += a.GetACCAgreementFromChain(c, dssp_assigned=True)
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
