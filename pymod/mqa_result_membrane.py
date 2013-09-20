from qmean import *
from qmean import mqa_membrane
from qmean import score_calculator
from qmean import reference_set
from qmean import conf
from ost.table import *
from ost.bindings import dssp
import os


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

  interaction   Normalized all-atom potential energy of the residue calculated from the
                short-range statistical potentials.

  cbeta         Normalized cbeta potential energy of the residue calculated from the
                short-range statistical potentials.

  packing       Normalized solvation potential energy of the residue calculated from
                the short-range statistical potentials.

  torsion       Torsion energy of the residue

  exposed       Relative solvent accessibility, calculated by dividing the maximally
                accessible surface area (ASA) of a residue by the observed value.
                ASA values are calculated by the DSSP program.

  ss_agreement  Agreement score of secondary structure predicted by psipred, with
                the observed seccondary structure observed in the model (DSSP).

  acc_agreement Agreement score of burial status predicted by accpro, with the
                the observed burial status in the structure (DSSP).

  lDDT          Predicted local lDDT value 
'''



GSCORES_TABLE_HEADER='''\
This file contains global scores calculated by the QMEAN scoring function.

References:

  Benkert, P., Tosatto, S.C.E. and Schomburg, D. (2008). "QMEAN: A comprehensive
  scoring function for model quality assessment." Proteins: Structure, Function,
  and Bioinformatics, 71(1):261-277.


  Benkert, P., Biasini, M. and Schwede, T. (2011). "Toward the estimation of the
  absolute quality of individual protein structure models." Bioinformatics
  doi: 10.1093/bioinformatics/btq662
'''


class MQAScore:
  def __init__(self, name, norm, z_score):
    self.name = name
    self.norm = norm
    self.z_score = z_score

  def __str__(self):
    return '%10s  norm=%.3f,  z_score=%.3f' % (self.name, 
                                               self.norm, 
                                               self.z_score)
  def __repr__(self):
    return 'MQAScore("%s",%f,%f)' % (self.name, self.norm, self.z_score)

class GlobalMembraneResult:
  def __init__(self, model, query, tab, ref_set):

    self.model = model
    self.ref_set = ref_set
    self.query = query
     
    try:
      name_idx = tab.GetColIndex('name')
      norm_idx = tab.GetColIndex('norm')
      z_idx = tab.GetColIndex('z_score')
    except:
      raise ValueError('The provided table must contain the columns name, norm and zscore!')
   
    required_features = ['QMEAN4','interaction','cbeta','packing','torsion']

    for f in required_features:
      if f not in [r[name_idx] for r in tab.rows]:
        raise ValueError('Provided table does not contain required feature %s!' % (f))

    for r in tab.rows:
      if r[name_idx] == 'QMEAN4':
        self.qmean4 = MQAScore('QMEAN4',r[norm_idx],r[z_idx])
      elif r[name_idx] == 'interaction':
        self.interaction = MQAScore('interaction',r[norm_idx],r[z_idx])
      elif r[name_idx] == 'cbeta':
        self.cbeta = MQAScore('cbeta',r[norm_idx],r[z_idx])
      elif r[name_idx] == 'packing':
        self.packing = MQAScore('packing',r[norm_idx],r[z_idx])
      elif r[name_idx] == 'torsion':
        self.torsion = MQAScore('torsion',r[norm_idx],r[z_idx])
      else:
        raise ValueError('Provided tab contains invalid feature %s!' % (r[name_idx]))


  def ReferenceSetPlot(self,feature='QMEAN4'):

    if feature == 'QMEAN4':
      return self.ref_set.ReferenceSetPlot(self.model.residue_count, self.qmean4.norm,'QMEAN4')
    else:
      raise ValueError('Can only create reference plots for QMEAN4!')

  def ZScoreSlider(self, score):
    return self.ref_set.ZScoreSliders([score.norm],[score.name],
                                       self.model.residue_count,compact=True)

  def ZScoreSliders(self, scores):
    norm_scores = list()
    names = list()
    labels = list()

    conversions = {'interaction':'Interaction', 'cbeta':'CBeta','packing':'Packing',
                   'torsion':'Torsion'}

    for s in scores:
      norm_scores.append(s.norm)
      names.append(s.name)
      if s.name in conversions:
        labels.append(conversions[s.name])
      else:
        labels.append(s.name)
    return self.ref_set.ZScoreSliders(norm_scores, names, self.model.residue_count,score_labels=labels)

  @staticmethod
  def Create(model, query, settings, psipred=None):

    tab=Table(['name','norm','z_score'],'sff')

    ref_set = reference_set.ReferenceSet(Table.Load(settings.reference_tab,format='pickle'))
    pot_soluble = PotentialContainer.Load(settings.global_potentials_soluble)
    pot_membrane = PotentialContainer.Load(settings.global_potentials_membrane)
    scorer = score_calculator.GlobalScorer.Load(settings.global_scorer)
    global_mqa = mqa_membrane.MembraneScores(model, model, pot_soluble, pot_membrane, query, psipred=psipred, 
                                             assign_dssp=False)


    features = ['interaction','cbeta','packing','torsion']


    global_mqa.CalculateScores(features)
    data = global_mqa.GetAVGData(features)

    s = scorer.GetGlobalScore('soluble', data)
    z = ref_set.ZScoreFromNormScore('QMEAN4',model.residue_count,s)

    tab.AddRow({'name':'QMEAN4','norm':s,'z_score':z})

    z = ref_set.ZScoreFromNormScore('interaction',model.residue_count, data['interaction'])
    tab.AddRow({'name':'interaction','norm':data['interaction'],'z_score':z})

    z = ref_set.ZScoreFromNormScore('cbeta',model.residue_count, data['cbeta'])
    tab.AddRow({'name':'cbeta','norm':data['cbeta'],'z_score':z})

    z = ref_set.ZScoreFromNormScore('packing',model.residue_count, data['packing'])
    tab.AddRow({'name':'packing','norm':data['packing'],'z_score':z})

    z = ref_set.ZScoreFromNormScore('torsion',model.residue_count, data['torsion'])
    tab.AddRow({'name':'torsion','norm':data['torsion'],'z_score':z})

    tab.comment=GSCORES_TABLE_HEADER
    return GlobalMembraneResult(model, query, tab, ref_set)

  def __str__(self):
    return '\n'.join([str(score) for score in self.all_scores])

  @property
  def all_scores(self):
    return(self.qmean4,self.interaction,self.cbeta,
           self.packing,self.torsion)

  @property
  def score_table(self):
    tab=Table(['name','norm','z_score'],'sff')
    for score in self.all_scores:
      tab.AddRow([score.name, score.norm, score.z_score])
    tab.comment=GSCORES_TABLE_HEADER
    return tab



class LocalMembraneResult:
  def __init__(self,model,data):
    self.model = model
    self.score_table = data

  def PlotlDDTProfile(self, chain=None):

    if chain == None:
      p = self.score_table.Plot('rnum','lDDT', style='-r',
                                   x_title='Residue Number', 
                                   y_title='Predicted lDDT',
                                   title='Local Model Quality',
                                   y_range = [-0.1,1.1], linewidth=2.0)
      return self._PlotMembraneProfile(p)
    else:
      tab = self.score_table.Filter(chain=chain)
      p = tab.Plot('rnum','lDDT', style='-r',
                       x_title='Residue Number', 
                       y_title='Predicted lDDT',
                       title='Local Model Quality',
                       y_range = [-0.1,1.1], linewidth=2.0)

      return self._PlotMembraneProfile(p)

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
  def Create(model, query, settings, psipred=None, accpro=None):

    pot_membrane = PotentialContainer.Load(settings.local_potentials_membrane)
    pot_soluble = PotentialContainer.Load(settings.local_potentials_soluble)
    scorer_membrane = score_calculator.LocalScorer.Load(settings.local_scorer_membrane)
    scorer_soluble = score_calculator.LocalScorer.Load(settings.local_scorer_soluble)
    local_mqa = mqa_membrane.MembraneScores(model, model, pot_soluble, pot_membrane, query, smooth_std=5.0, psipred=psipred, 
                            accpro=accpro, assign_dssp=False)

    local_mqa.CalculateScores(['interaction','cbeta','packing','torsion','exposed'])
    data = local_mqa.GetLocalData(['interaction','cbeta','packing','torsion','exposed'])

    for r in model.residues:
      r.SetIntProp('in_membrane',0)

    membrane_view = model.Select(query)
    for r in membrane_view.residues:
      r.SetIntProp('in_membrane',1)


    try:
      local_mqa.CalculateScores(['ss_agreement'])
      temp = local_mqa.GetLocalData(['ss_agreement'])['ss_agreement']
      for i,r in enumerate(model.residues):
        if r.GetIntProp('in_membrane')==1:
          temp[i] = float('NaN')
      data['ss_agreement'] = temp
    except:
      data['ss_agreement'] = model.residue_count*[float('NaN')]
    try:
      local_mqa.CalculateScores(['acc_agreement'])
      temp = local_mqa.GetLocalData(['acc_agreement'])['acc_agreement']
      for i,r in enumerate(model.residues):
        if r.GetIntProp('in_membrane')==1:
          temp[i] = float('NaN')
      data['acc_agreement'] = temp['acc_agreement']
    except:
      data['acc_agreement'] = model.residue_count*[float('NaN')]


    scores = list()
    dssp_ss = local_mqa.GetDSSPSS()

    for i ,r in enumerate(model.residues):
      if r.GetIntProp('in_membrane')==0:
        residue_data=dict()
        for f in ['interaction','cbeta','packing','torsion','exposed','ss_agreement','acc_agreement']:
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
        for f in ['interaction','cbeta','packing','torsion','exposed']:
          residue_data[f] = data[f][i]
        if local_mqa.transmembrane_type=='helical':
          scores.append(scorer_membrane.GetLocalScore('helical', r.one_letter_code, residue_data))
        else:
          scores.append(scorer_membrane.GetLocalScore('extended', r.one_letter_code, residue_data))

    data['lDDT'] = scores

    lscores=Table(['chain', 'rindex', 'rnum', 'rname', 'membrane','interaction',
                   'cbeta', 'packing', 'torsion', 'exposed',
                   'ss_agreement', 'acc_agreement', 'lDDT'],
                   'siisiffffffff')

    for i, res in enumerate(model.residues):
      lscores.AddRow({ 'chain' : res.chain.name,
                      'rindex' : res.index,
                      'rnum' : res.number.num,
                      'rname' : res.name,
                      'membrane' : res.GetIntProp('in_membrane'),
                      'interaction' : data['interaction'][i],
                      'cbeta' : data['cbeta'][i],
                      'packing' : data['packing'][i],
                      'torsion' : data['torsion'][i],
                      'exposed' : data['exposed'][i],
                      'ss_agreement' : data['ss_agreement'][i],
                      'acc_agreement' : data['acc_agreement'][i],
                      'lDDT' : data['lDDT'][i]})

    lscores.comment=LSCORES_TABLE_HEADER
    return LocalMembraneResult(model,lscores)



def AssessModelQuality(model, query, output_dir='.', plots=True, local_scores=True,
                       global_scores=True, table_format='ost', psipred=None, 
                       accpro=None,settings=conf.MembraneSettings(), dssp_path=None):
  #hack, that torsion handles get assigned
  processor = conop.HeuristicProcessor(connect=False,peptide_bonds=False,assign_torsions=True)
  processor.Process(model.handle,False)

  dssp.AssignDSSP(model, extract_burial_status=True,dssp_bin=dssp_path)

  results = []

  if not os.path.exists(output_dir):
    os.makedirs(output_dir)
  plot_dir=os.path.join(output_dir,'plots')
  sliders_dir=os.path.join(plot_dir,'sliders')
  if plots:
    if not os.path.exists(sliders_dir):
      os.makedirs(sliders_dir)
  if global_scores:
    global_result=GlobalMembraneResult.Create(model,query,settings, psipred=psipred)
    tab=global_result.score_table
    tab.Save(os.path.join(output_dir,'global_scores.txt'),format=table_format)
    if plots:
      p=global_result.ReferenceSetPlot('QMEAN4')
      p.savefig(os.path.join(plot_dir,'qmean4_reference_set.png'))
      qmean4_scores = [global_result.qmean4, global_result.interaction, global_result.cbeta, global_result.packing, global_result.torsion]
      p=global_result.ZScoreSliders(qmean4_scores)
      p.savefig(os.path.join(plot_dir,'qmean4_slider.png'))
      for s in qmean4_scores:
        p=global_result.ZScoreSlider(s)
        p.savefig(os.path.join(sliders_dir,'%s_compact.png'%(s.name)))
    results.append(global_result)
  if local_scores:
    local_result=LocalMembraneResult.Create(model,query,settings,psipred=psipred,accpro=accpro)
    tab=local_result.score_table
    tab.Save(os.path.join(output_dir,'local_scores.txt'),format=table_format)
    if plots:
      for ch in model.chains:
        p=local_result.PlotlDDTProfile(chain=ch.name)
        p.savefig(os.path.join(plot_dir,'local_lDDT_%s.png' % (ch.name)))
    results.append(local_result)

  return results





