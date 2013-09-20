from qmean import *
from qmean import conf
from ost.bindings import dssp
from ost import mol


class MembraneScores:

  def __init__(self, target, environment, potential_container_soluble, potential_container_membrane, membrane_query, smooth_std=None, psipred=None, accpro=None, assign_dssp=True, norm=True):

    self.data=dict()
    self.target=target
    self.environment=environment
    self.potential_container_soluble=potential_container_soluble
    self.potential_container_membrane=potential_container_membrane
    self.smooth_std=smooth_std
    self.psipred=psipred
    self.accpro=accpro
    self.norm=norm
    self.membrane_query=membrane_query
    self.ca_positions=None
    self.spherical_smoother=None
    self.transmembrane_type=None

    self.membrane_target=target.Select(membrane_query)
    self.soluble_target=mol.Difference(target, self.membrane_target)
    self.membrane_environment=environment.Select(membrane_query)
    self.soluble_environment=mol.Difference(environment, self.membrane_environment)
    
    if target.handle != environment.handle:
      raise RuntimeError('Target and Environment must be views of the same handle!')

    if assign_dssp:
      dssp.AssignDSSP(target, extract_burial_status=True)      

    #if there is a membrane, check whether the structure has a helical or extended
    #transmembranepart

    num_helical_membrane_residues=0
    num_extended_membrane_residues=0

    for r in self.membrane_target.residues:
      if r.GetSecStructure().IsHelical():
        num_helical_membrane_residues+=1
      elif r.GetSecStructure().IsExtended():
        num_extended_membrane_residues+=1

    self.transmembrane_type=""
    if num_helical_membrane_residues > num_extended_membrane_residues:
      self.transmembrane_type="helical"
    else:
      self.transmembrane_type="extended"


    #if psipred data is provided, take secondary information from there. Particularly in models with
    #low quality, the secondary structure assignment of DSSP srews up. In case of membrane proteins,
    #the secondary structure of the membrane spanning region is either set to helical or extended, 
    #depending on the majority of residues.

    self.dssp_ss = list()

    for r in target.residues:
      if r.GetSecStructure().IsHelical():
        self.dssp_ss.append('H')
      if r.GetSecStructure().IsExtended():
        self.dssp_ss.append('E')
      else:
        self.dssp_ss.append('C')


    if psipred != None:
      temp_ss = psipred.GetPSIPREDSS(target)
    else:
      temp_ss = self.dssp_ss


    self.sec_structure=list()
    for ss, r in zip(temp_ss, self.target.residues):
      if self.membrane_target.ViewForHandle(r.handle).IsValid():
        if self.transmembrane_type=='helical':
          self.sec_structure.append('H')
        else:
          self.sec_structure.append('E')
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
    self.soluble_helix_selection = self.soluble_target.Select('grSS=0')
    self.soluble_extended_selection = self.soluble_target.Select('grSS=1')
    self.soluble_coil_selection = self.soluble_target.Select('grSS=2')

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

  def GetSolubleTarget(self):
    return self.soluble_target

  def GetMembraneTarget(self):
    return self.membrane_target

  def GetSolubleEnvironment(self):
    return self.soluble_environment

  def GetMembraneEnvironment(self):
    return self.membrane_environment

  def GetDSSPSS(self):
    return self.dssp_ss

  def GetAssignedSS(self):
    return self.sec_structure

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
          raise ValueError("you must provide accpro data for acc_agreement term!")
      elif f == 'fraction_loops':
        self.GetFractionLoops()
      elif f == 'exposed':
        self.GetExposed()
      else:
        raise ValueError('Requested feature is not supported in local scorer!')

  
  def GetTorsion(self):


    if self.smooth_std==None:
      soluble_torsion = self.potential_container_soluble['torsion'].GetEnergies(self.soluble_target, normalize=self.norm)
      if self.transmembrane_type=="helical":
        membrane_torsion = self.potential_container_membrane['torsion_helix'].GetEnergies(self.membrane_target, normalize=self.norm)
      else:
        membrane_torsion = self.potential_container_membrane['torsion_extended'].GetEnergies(self.membrane_target, normalize=self.norm)

      torsion_energies=list()
      for r in self.target.residues:
        if self.membrane_target.ViewForHandle(r.handle).IsValid():
          torsion_energies.append(membrane_torsion.pop(0))
        else:
          torsion_energies.append(soluble_torsion.pop(0))
      self.data['torsion'] = torsion_energies

      overall_soluble_torsion = self.potential_container_soluble['torsion'].GetEnergies(self.target,normalize=self.norm)
      self.data['avg_torsion'] = self.GetAverage(overall_soluble_torsion)

    else:
      soluble_torsion_helical=self.potential_container_soluble['simple_torsion_helix'].GetEnergies(self.soluble_helix_selection, normalize=self.norm)
      soluble_torsion_extended=self.potential_container_soluble['simple_torsion_extended'].GetEnergies(self.soluble_extended_selection, normalize=self.norm)
      soluble_torsion_coil=self.potential_container_soluble['simple_torsion_coil'].GetEnergies(self.soluble_coil_selection, normalize=self.norm)
      if self.transmembrane_type=="helical":
        membrane_torsion=self.potential_container_membrane['torsion_helix'].GetEnergies(self.membrane_target, normalize=self.norm)
      else:
        membrane_torsion=self.potential_container_membrane['torsion_extended'].GetEnergies(self.membrane_target, normalize=self.norm)


      temp_membrane=list()

      for r in self.target.residues:
        if self.membrane_target.ViewForHandle(r.handle).IsValid():
          temp_membrane.append(membrane_torsion.pop(0))
        else:
          if r.GetIntProp('SS')==0:
            temp_membrane.append(soluble_torsion_helical.pop(0))
          elif r.GetIntProp('SS')==1:
            temp_membrane.append(soluble_torsion_extended.pop(0))
          else:
            temp_membrane.append(soluble_torsion_coil.pop(0))

      smoothed_temp_membrane = self.spherical_smoother.Smooth(temp_membrane)

      temp_soluble = self.potential_container_soluble['torsion'].GetEnergies(self.target,normalize=self.norm)
      smoothed_temp_soluble = self.spherical_smoother.Smooth(temp_soluble)

      torsion_energies=list()
      for r, s_t_m, s_t_s in zip(self.target.residues, smoothed_temp_membrane, smoothed_temp_soluble):
        if self.membrane_target.ViewForHandle(r.handle).IsValid():
          torsion_energies.append(s_t_m)
        else:
          torsion_energies.append(s_t_s)

      self.data['torsion'] = torsion_energies

      self.data['avg_torsion'] = self.GetAverage(temp_soluble)
 
  def GetInteraction(self):


    interaction=list()

    soluble_helix_interaction = self.potential_container_soluble['interaction_helix'].GetEnergies(self.soluble_helix_selection, self.soluble_environment, normalize=self.norm)
    soluble_extended_interaction = self.potential_container_soluble['interaction_extended'].GetEnergies(self.soluble_extended_selection, self.soluble_environment, normalize=self.norm)
    soluble_coil_interaction = self.potential_container_soluble['interaction_coil'].GetEnergies(self.soluble_coil_selection, self.soluble_environment, normalize=self.norm)

    if self.transmembrane_type=="helical":
      membrane_interaction = self.potential_container_membrane['interaction_helix'].GetEnergies(self.membrane_target,self.membrane_environment,normalize=self.norm)
    else:
      membrane_interaction = self.potential_container_membrane['interaction_extended'].GetEnergies(self.membrane_target,self.membrane_environment, normalize=self.norm)

    for r in self.target.residues:
      if self.membrane_target.ViewForHandle(r.handle).IsValid():
        interaction.append(membrane_interaction.pop(0))
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

    soluble_helix_cbeta = self.potential_container_soluble['cbeta_helix'].GetEnergies(self.soluble_helix_selection, self.soluble_environment, normalize=self.norm)
    soluble_extended_cbeta = self.potential_container_soluble['cbeta_extended'].GetEnergies(self.soluble_extended_selection, self.soluble_environment, normalize=self.norm)
    soluble_coil_cbeta = self.potential_container_soluble['cbeta_coil'].GetEnergies(self.soluble_coil_selection, self.soluble_environment, normalize=self.norm)

    if self.transmembrane_type=="helical":
      membrane_cbeta = self.potential_container_membrane['cbeta_helix'].GetEnergies(self.membrane_target, self.membrane_environment, normalize=self.norm)
    else:
      membrane_cbeta = self.potential_container_membrane['cbeta_extended'].GetEnergies(self.membrane_target, self.membrane_environment, normalize=self.norm)

    for r in self.target.residues:
      if self.membrane_target.ViewForHandle(r.handle).IsValid():
        cbeta.append(membrane_cbeta.pop(0))
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


  def GetReduced(self):
    reduced=list()

    soluble_helix_reduced = self.potential_container_soluble['reduced_helix'].GetEnergies(self.soluble_helix_selection, self.soluble_environment, normalize=self.norm)
    soluble_extended_reduced = self.potential_container_soluble['reduced_extended'].GetEnergies(self.soluble_extended_selection, self.soluble_environment, normalize=self.norm)
    soluble_coil_reduced = self.potential_container_soluble['reduced_coil'].GetEnergies(self.soluble_coil_selection, self.soluble_environment, normalize=self.norm)

    if self.transmembrane_type=="helical":
      membrane_reduced = self.potential_container_membrane['reduced_helix'].GetEnergies(self.membrane_target, self.membrane_environment, normalize=self.norm)
    else:
      membrane_reduced = self.potential_container_membrane['reduced_extended'].GetEnergies(self.membrane_target, self.membrane_environment, normalize=self.norm)

    for r in self.target.residues:
      if self.membrane_target.ViewForHandle(r.handle).IsValid():
        reduced.append(membrane_reduced.pop(0))
      else:
        if r.GetIntProp('SS')==0:
          reduced.append(soluble_helix_reduced.pop(0))
        elif r.GetIntProp('SS')==1:
          reduced.append(soluble_extended_reduced.pop(0))
        else:
          reduced.append(soluble_coil_reduced.pop(0))

    if self.smooth_std!=None:
      self.data['reduced'] = self.spherical_smoother.Smooth(reduced)
    else:
      self.data['reduced'] = reduced

    self.data['avg_reduced'] = self.GetAverage(reduced)


  def GetPacking(self):
    packing=list()

    soluble_helix_packing = self.potential_container_soluble['packing_helix'].GetEnergies(self.soluble_helix_selection, self.environment,normalize=self.norm)
    soluble_extended_packing = self.potential_container_soluble['packing_extended'].GetEnergies(self.soluble_extended_selection, self.environment,normalize=self.norm)
    soluble_coil_packing = self.potential_container_soluble['packing_coil'].GetEnergies(self.soluble_coil_selection, self.environment,normalize=self.norm)

    if self.transmembrane_type=="helical":
      membrane_packing = self.potential_container_membrane['packing_helix'].GetEnergies(self.membrane_target, self.environment, normalize=self.norm)
    else:
      membrane_packing = self.potential_container_membrane['packing_extended'].GetEnergies(self.membrane_target, self.environment, normalize=self.norm)

    for r in self.target.residues:
      if self.membrane_target.ViewForHandle(r.handle).IsValid():
        packing.append(membrane_packing.pop(0))
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
        exposed.append(r.GetFloatProp('relative_solvent_accessibility'))
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
      ss_agreement=self.psipred.GetSSAgreementFromChain(self.full_target, dssp_assigned=True)
      if self.smooth_std!=None:
        self.data['ss_agreement']=self.spherical_smoother.Smooth(ss_agreement)
      else:
        self.data['ss_agreement']=ss_agreement
      self.data['avg_ss_agreement'] = self.GetAverage(ss_agreement)
    else:
      raise ValueError("Cannot calculate ss_agreement term without psipred information!")

  def GetACCAgreement(self):
    if self.accpro!=None:
      acc_agreement=self.accpro.GetACCAgreementFromChain(self.full_target, dssp_assigned=True)
      if self.smooth_std!=None:
        self.data['acc_agreement']=self.spherical_smoother.Smooth(acc_agreement)
      else:
        self.data['acc_agreement'] = acc_agreement
      self.data['avg_acc_agreement'] = self.GetAverage(acc_agreement)
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
