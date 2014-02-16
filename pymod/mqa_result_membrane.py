from qmean import *
from qmean import mqa_membrane_hybrid
from qmean import score_calculator
from qmean import reference_set
from qmean import conf
from ost.table import *
from ost.bindings import dssp
import os
import matplotlib.pyplot as plt


LSCORES_TABLE_HEADER='''\
This file contains local scores calculated by the QMEAN scoring function.

References:

  Benkert, P., Tosatto, S.C.E. and Schomburg, D. (2008). "QMEAN: A comprehensive
  scoring function for model quality assessment." Proteins: Structure, Function,
  and Bioinformatics, 71(1):261-277.


  Benkert, P., Biasini, M. and Schwede, T. (2011). "Toward the estimation of the
  absolute quality of individual protein structure models." Bioinformatics
  doi: 10.1093/bioinformatics/btq662

Description of columns:

  chain         Chainname as written in the PDB file
  rindex        Zero-base index of residue in the chain
  rnum          Residue number as written in the PDB file
  rname         Residue three-letter-code, e.g. ARG

  all_atom      Normalized all-atom potential energy of the residue calculated from the
                short-range statistical potentials.

  cbeta         Normalized cbeta potential energy of the residue calculated from the
                short-range statistical potentials.

  solvation     Normalized solvation potential energy of the residue calculated from
                the short-range statistical potentials.

  torsion       Torsion energy of the residue

  exposed       Relative solvent accessibility, calculated by dividing the maximally
                accessible surface area (ASA) of a residue by the observed value.
                ASA values are calculated by the DSSP program.

  ss_agreement  Agreement score of secondary structure predicted by psipred, with
                the observed seccondary structure observed in the model (DSSP).

  acc_agreement Agreement score of burial status predicted by accpro, with the
                the observed burial status in the structure (DSSP).

  membrane      Assignment of the membrane part. 1 for membrane residues, 0 otherwise

  lDDT          Predicted local lDDT value 
'''


class LocalMembraneResult:
  def __init__(self,model,data):
    self.model = model
    self.score_table = data

  def PlotlDDTProfile(self, chain=None):

    chain_score_table = self.score_table
    if chain != None and chain !='':
      chain_score_table = self.score_table.Filter(chain=chain)

    rnums = list()
    scores = list()
    rnum_idx = chain_score_table.GetColIndex('rnum')
    score_idx = chain_score_table.GetColIndex('lDDT')
    for r in chain_score_table.rows:
      rnums.append(r[rnum_idx])
      scores.append(r[score_idx])

    rnum_slices = list()
    score_slices = list()
    current_rnum_slice = list()
    current_score_slice = list()

    previous_rnum = rnums[0]-1

    for i in range(len(rnums)):
      if rnums[i] != previous_rnum + 1:
        rnum_slices.append(current_rnum_slice)
        score_slices.append(current_score_slice)
        current_rnum_slice = list()
        current_score_slice = list()
      current_rnum_slice.append(rnums[i])
      current_score_slice.append(scores[i])
      previous_rnum = rnums[i]

    score_slices.append(current_score_slice)
    rnum_slices.append(current_rnum_slice)

    plt.clf()
    for s,r in zip(score_slices, rnum_slices):
      plt.plot(r,s,color='r',linewidth=2.0)

    plt.ylim([0,1])
    plt.xlim([min(rnums),max(rnums)])
    
    return plt

  def _PlotMembraneProfile(self, p):
    
    membrane_idx = self.score_table.GetColIndex('membrane')
    rnum_idx = self.score_table.GetColIndex('rnum')
    starts = list()
    ends = list()
    in_membrane = False

    old_row = None

    for r in self.score_table.rows:
      if r[membrane_idx] not in [0,1]:
        raise ValueError('Invalid Membrane Association Value!')
      if in_membrane==False and r[membrane_idx]==1:
        starts.append(r[rnum_idx])
        in_membrane=True
      elif in_membrane==True and r[membrane_idx]==0:
        ends.append(old_row[rnum_idx])
        in_membrane=False
      old_row=r

    if len(starts)!=len(ends):
      ends.append(self.score_table.rows[-1][rnum_idx])

    if len(starts)!=len(ends):
      raise ValueError('Error Plotting Local Membrane Values!')


    print starts,ends

    for s,e in zip(starts, ends):
      p.axvspan(s, e, color='#004586',alpha=0.2)

    return p

  @staticmethod
  def Create(model, settings, membrane_query = None, interface_query = None, mem_param = None, psipred=None, accpro=None):

    pot_membrane = PotentialContainer.Load(settings.local_potentials_membrane)
    pot_soluble = PotentialContainer.Load(settings.local_potentials_soluble)
    scorer_membrane = score_calculator.LocalScorer.Load(settings.local_scorer_membrane)
    scorer_soluble = score_calculator.LocalScorer.Load(settings.local_scorer_soluble)
    local_mqa = mqa_membrane_hybrid.MembraneScores(model, model, pot_soluble, pot_membrane, smooth_std=5.0, psipred=psipred, 
                            accpro=accpro, membrane_query = membrane_query, interface_query = interface_query,mem_param = mem_param)


    features = ['interaction','cbeta','packing','torsion','exposed']
    if psipred != None:
      features.append('ss_agreement')
    if accpro != None:
      features.apend('acc_agreement')

    local_mqa.CalculateScores(features)
    data = local_mqa.GetLocalData(features)

    membrane_states = local_mqa.GetMembraneState()

    scores = list()
    dssp_ss = local_mqa.GetDSSPSS()

    for i ,r in enumerate(model.residues):
      if membrane_states[i] == 0:
        residue_data=dict()
        for f in features:
          residue_data[f]=data[f][i]
        ss = ''
        if dssp_ss[i]=='H':
          ss = 'helical'
        elif dssp_ss[i]=='E':
          ss = 'extended'
        else:
          ss = 'coil'
        scores.append(scorer_soluble.GetLocalScore(ss, r.one_letter_code, residue_data))

      else:
        residue_data=dict()
        for f in features:
          if f in ['ss_agreement','acc_agreement']:
            continue
          residue_data[f] = data[f][i]
        if membrane_states[i] == 1:
          scores.append(scorer_membrane.GetLocalScore('membrane', r.one_letter_code, residue_data))
        elif membrane_states[i] == 2:
          scores.append(scorer_membrane.GetLocalScore('interface', r.one_letter_code, residue_data))
        else:
          raise RuntimeError("Invalid membrane state encountered")

    data['lDDT'] = scores
    #check whether ss_agreement and acc_agreement are set, otherwise we just
    #add NaN values
    if 'ss_agreement' not in features:
      data['ss_agreement'] = [float('NaN')] * len(data[features[0]])
    if 'acc_agreement' not in features:
      data['acc_agreement'] = [float('NaN')] * len(data[features[0]])

    lscores=Table(['chain', 'rindex', 'rnum', 'rname', 'membrane','all_atom',
                   'cbeta', 'solvation', 'torsion', 'exposed',
                   'ss_agreement', 'acc_agreement', 'lDDT'],
                   'siisiffffffff')

    for i, res in enumerate(model.residues):
      lscores.AddRow({ 'chain' : res.chain.name,
                      'rindex' : res.index,
                      'rnum' : res.number.num,
                      'rname' : res.name,
                      'membrane' : res.GetIntProp('in_membrane'),
                      'all_atom' : data['interaction'][i],
                      'cbeta' : data['cbeta'][i],
                      'solvation' : data['packing'][i],
                      'torsion' : data['torsion'][i],
                      'exposed' : data['exposed'][i],
                      'ss_agreement' : data['ss_agreement'][i],
                      'acc_agreement' : data['acc_agreement'][i],
                      'lDDT' : data['lDDT'][i]})

    lscores.comment=LSCORES_TABLE_HEADER
    return LocalMembraneResult(model,lscores)


def AssessModelQuality(model, membrane_query=None, interface_query=None, mem_param = None, output_dir='.', plots=True,
                       table_format='ost', psipred=None, accpro=None,settings=conf.MembraneSettings(), 
                       dssp_path=None):
  #hack, that torsion handles get assigned
  processor = conop.HeuristicProcessor(connect=False,peptide_bonds=False,assign_torsions=True)
  processor.Process(model.handle,False)

  dssp.AssignDSSP(model, extract_burial_status=True,dssp_bin=dssp_path)

  if not os.path.exists(output_dir):
    os.makedirs(output_dir)
  plot_dir=os.path.join(output_dir,'plots')
  if plots:
    if not os.path.exists(plot_dir):
      os.makedirs(plot_dir)
  local_result=LocalMembraneResult.Create(model,settings,membrane_query,interface_query,mem_param,psipred=psipred,accpro=accpro)
  tab=local_result.score_table
  tab.Save(os.path.join(output_dir,'local_scores.txt'),format=table_format)
  if plots:
    for ch in model.chains:
      p=local_result.PlotlDDTProfile(chain=ch.name)
      p.savefig(os.path.join(plot_dir,'local_lDDT_%s.png' % (ch.name)))
  return local_result





