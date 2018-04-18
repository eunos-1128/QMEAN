from qmean import *
from qmean import mqa
from qmean import score_calculator
from qmean import reference_set
from qmean import conf
from ost.table import *
from ost.bindings import dssp
from sklearn.ensemble import RandomForestRegressor
import os, pickle
import numpy as np


sm_conversions = {'QMEAN4':'qmean4','QMEAN6':'qmean6','interaction':'all_atom','cbeta':'cbeta','packing':'solvation','torsion':'torsion','ss_agreement':'ss_agreement','acc_agreement':'acc_agreement'}
label_conversions = {'interaction':'All Atom', 'cbeta':'CBeta','packing':'Solvation',
                     'torsion':'Torsion','ss_agreement':'SS Agree','acc_agreement':'ACC Agree',
                     'QMEAN4':'QMEAN4','QMEAN6':'QMEAN6'}

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

  dist_const    Score predicting correctness of pairwise residue distances in the model.
                Distance constraints are obtained from distances between residue pairs
                in homologous templates.

  QMEAN         Predicted local quality. A value describing the expected local
                similarity to the target structure with a range of [0,1]. 

  QMEANDisCo    Predicted local quality by combining QMEAN, dist_const and sequence 
                based features using a random forest. The range remains [0,1].  
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


class GlobalResult:
  def __init__(self, model, tab, ref_set, global_mqa):

    self.model = model
    self.ref_set = ref_set
    self.global_mqa = global_mqa
     
    try:
      name_idx = tab.GetColIndex('name')
      norm_idx = tab.GetColIndex('norm')
      z_idx = tab.GetColIndex('z_score')
    except:
      raise ValueError('The provided table must contain the columns name, norm and zscore!')
   
    required_features = ['qmean4','qmean6','all_atom','cbeta','solvation','torsion','ss_agreement','acc_agreement']

    for f in required_features:
      if f not in [r[name_idx] for r in tab.rows]:
        raise ValueError('Provided table does not contain required feature %s!' % (f))

    for r in tab.rows:
      if r[name_idx] == 'qmean4':
        self.qmean4 = MQAScore('QMEAN4',r[norm_idx],r[z_idx])
      elif r[name_idx] == 'qmean6':
        self.qmean6 = MQAScore('QMEAN6',r[norm_idx],r[z_idx])
      elif r[name_idx] == 'all_atom':
        self.interaction = MQAScore('interaction',r[norm_idx],r[z_idx])
      elif r[name_idx] == 'cbeta':
        self.cbeta = MQAScore('cbeta',r[norm_idx],r[z_idx])
      elif r[name_idx] == 'solvation':
        self.packing = MQAScore('packing',r[norm_idx],r[z_idx])
      elif r[name_idx] == 'torsion':
        self.torsion = MQAScore('torsion',r[norm_idx],r[z_idx])
      elif r[name_idx] == 'ss_agreement':
        self.ss_agreement = MQAScore('ss_agreement',r[norm_idx],r[z_idx])
      elif r[name_idx] == 'acc_agreement':
        self.acc_agreement = MQAScore('acc_agreement',r[norm_idx],r[z_idx])
      else:
        raise ValueError('Provided tab contains invalid feature %s!' % (r[name_idx]))


  def ReferenceSetPlot(self,feature='qmean4'):

    if feature.lower() == 'qmean4':
      return self.ref_set.ReferenceSetPlot(self.model.residue_count, self.qmean4.norm,'QMEAN4')
    elif feature.lower() =='qmean6':
      return self.ref_set.ReferenceSetPlot(self.model.residue_count, self.qmean6.norm,'QMEAN6')
    else:
      raise ValueError('Can only create reference plots for qmean4 or qmean6!')

  def ZScoreSlider(self, score, compact=True):
    label = [label_conversions[score.name]]
    return self.ref_set.ZScoreSliders([score.norm],[score.name], self.model.residue_count,
                                      compact=compact,score_labels = label)

  def ZScoreSliders(self, scores):
    norm_scores = list()
    names = list()
    labels = list()

    for s in scores:
      norm_scores.append(s.norm)
      names.append(s.name)
      if s.name in label_conversions:
        labels.append(label_conversions[s.name])
      else:
        labels.append(s.name)
    return self.ref_set.ZScoreSliders(norm_scores, names, self.model.residue_count,score_labels=labels)

  @staticmethod
  def Create(model, settings, psipred=None, accpro=None):

    tab=Table(['name','norm','z_score'],'sff')

    ref_set = reference_set.ReferenceSet(Table.Load(settings.reference_tab))
    pot = PotentialContainer.Load(settings.global_potentials)
    scorer = score_calculator.GlobalScorer.Load(settings.global_scorer)
    global_mqa = mqa.Scores(model, model, pot, psipred=psipred, 
                            accpro=accpro, assign_dssp=False)

    agreement_terms = False

    if psipred!=None and accpro!=None:
      agreement_terms=True

    if agreement_terms:
      features = ['interaction','cbeta','packing','torsion','ss_agreement','acc_agreement']
    else:
      features = ['interaction','cbeta','packing','torsion']

    data = global_mqa.GetAVGData(['interaction','cbeta','torsion','packing'])
    s = min(max(scorer.GetGlobalScore('soluble', data),0.0),1.0)
    z = ref_set.ZScoreFromNormScore('QMEAN4',model.residue_count,s)

    tab.AddRow({'name':'qmean4','norm':s,'z_score':z})
    
    if agreement_terms:
      data = global_mqa.GetAVGData(['interaction','cbeta','torsion','packing','ss_agreement','acc_agreement'])
      s = scorer.GetGlobalScore('soluble',data)
      z = ref_set.ZScoreFromNormScore('QMEAN6',model.residue_count,s)
      tab.AddRow({'name':'qmean6','norm':s,'z_score':z})
    else:
      tab.AddRow({'name':'qmean6', 'norm':float('NaN'), 'z_score':float('NaN')})

    z = ref_set.ZScoreFromNormScore('interaction',model.residue_count, data['interaction'])
    tab.AddRow({'name':'all_atom','norm':data['interaction'],'z_score':z})

    z = ref_set.ZScoreFromNormScore('cbeta',model.residue_count, data['cbeta'])
    tab.AddRow({'name':'cbeta','norm':data['cbeta'],'z_score':z})

    z = ref_set.ZScoreFromNormScore('packing',model.residue_count, data['packing'])
    tab.AddRow({'name':'solvation','norm':data['packing'],'z_score':z})

    z = ref_set.ZScoreFromNormScore('torsion',model.residue_count, data['torsion'])
    tab.AddRow({'name':'torsion','norm':data['torsion'],'z_score':z})

    if agreement_terms:
      z = ref_set.ZScoreFromNormScore('ss_agreement',model.residue_count, data['ss_agreement'])
      tab.AddRow({'name':'ss_agreement','norm':data['ss_agreement'],'z_score':z})

      z = ref_set.ZScoreFromNormScore('acc_agreement',model.residue_count, data['acc_agreement'])
      tab.AddRow({'name':'acc_agreement','norm':data['acc_agreement'],'z_score':z})
    else:
      tab.AddRow({'name':'ss_agreement','norm':float('NaN'),'z_score':float('NaN')})
      tab.AddRow({'name':'acc_agreement','norm':float('NaN'),'z_score':float('NaN')})

    tab.comment=GSCORES_TABLE_HEADER
    return GlobalResult(model, tab, ref_set, global_mqa)

  def __str__(self):
    return '\n'.join([str(score) for score in self.all_scores])

  @property
  def all_scores(self):
    return(self.qmean4,self.qmean6,self.interaction,self.cbeta,
           self.packing,self.torsion,self.ss_agreement, self.acc_agreement)

  @property
  def score_table(self):
    tab=Table(['name','norm','z_score'],'sff')
    for score in self.all_scores:
      tab.AddRow([sm_conversions[score.name], score.norm, score.z_score])
    tab.comment=GSCORES_TABLE_HEADER
    return tab



class LocalResult:
  def __init__(self, model, data, local_mqa):
    self.model = model
    self.score_table = data
    self.local_mqa = local_mqa

  def PlotlDDTProfile(self, chain=None, plot_disco=False):

    from matplotlib import pyplot
    pyplot.clf()
    pyplot.figure(figsize=[8,6],dpi=80)
    pyplot.ylim((-0.1,1.1))
    pyplot.xlabel('Residue Number',size='large')
    pyplot.ylabel('Predicted Local Similarity to Target',size='large')

    color_scheme = ['#ff420e','#004586','#ffd320','#578d1c','#7e0021']

    chain_idx = self.score_table.GetColIndex('chain')
    rnum_idx = self.score_table.GetColIndex('rnum')
    qmean_idx = None
    if plot_disco:
      qmean_idx = self.score_table.GetColIndex('QMEANDisCo')
    else:
      qmean_idx = self.score_table.GetColIndex('QMEAN')


    chains = list()
    chains.append(self.score_table.rows[0][chain_idx])

    for r in self.score_table.rows:
      if r[chain_idx] != chains[-1]:
        chains.append(r[chain_idx])


    if chain == None:

      pyplot.title('Local Quality Estimate',size='x-large')

      for i,ch in enumerate(chains):
        chain_tab = self.score_table.Filter(chain=ch)  
        res_num = list()
        predicted_lddt = list()
        for r in chain_tab.rows:
          res_num.append(r[rnum_idx])
          predicted_lddt.append(r[qmean_idx])
        pyplot.plot(res_num,predicted_lddt,color=color_scheme[i%len(color_scheme)],linewidth=2.0)
      return pyplot
         
    else:

      pyplot.title('Local Quality Estimate: Chain %s'%(chain),size='x-large')
      chain_tab = self.score_table.Filter(chain=chain)
      res_num = list()
      predicted_lddt = list()
      for r in chain_tab.rows:
        res_num.append(r[rnum_idx])
        predicted_lddt.append(r[qmean_idx])
      color_idx = chains.index(chain)
      pyplot.plot(res_num,predicted_lddt,color=color_scheme[color_idx%len(color_scheme)],linewidth=2.0)

      return pyplot

  @staticmethod
  def Create(model, settings, assign_bfactors, 
             psipred = None, accpro = None, dc = None, global_result = None):

    pot = PotentialContainer.Load(settings.local_potentials)
    scorer = score_calculator.NNScorer.Load(settings.local_scorer)
    local_mqa = mqa.Scores(model, model, pot, smooth_std=5.0, psipred=psipred, 
                           accpro=accpro, dc=dc, assign_dssp=False)

    if global_result == None:
      global_result = GlobalResult.Create(model, settings, 
                                          psipred = psipred, accpro = accpro) 
    global_mqa = global_result.global_mqa

    global_features = ['torsion', 'reduced', 'interaction', 'cbeta', 
                       'packing', 'cb_packing']

    local_features = ['counts', 'packing', 'cb_packing', 'torsion', 'reduced', 
                      'interaction', 'cbeta', 'ss_agreement', 'acc_agreement', 
                      'dist_const']

    data = local_mqa.GetLocalData(local_features)
    global_data = global_mqa.GetAVGData(global_features)

    qmean_scores = list()
    for i, r in enumerate(model.residues):

      f_dict = dict()

      # do the global features
      for f in global_features:
        f_dict["avg_" + f] = global_data[f]

      # do the local features
      for f in local_features:
        # dist const requires special treatment
        if f == "dist_const":
          continue
        f_dict[f] = data[f][i]

      # do dist_const
      dist_const = data["dist_const"]
      s = dist_const["disco"][i]
      if s == s:
        f_dict["disco"] = s
        f_dict["disco_counts"] = dist_const["counts"][i]
        f_dict["disco_avg_num_clusters"] = dist_const["avg_num_clusters"][i]
        f_dict["disco_avg_max_seqsim"] = dist_const["avg_max_seqsim"][i]
        f_dict["disco_avg_max_seqid"] = dist_const["avg_max_seqid"][i]
        f_dict["disco_avg_variance"] = dist_const["avg_variance"][i]
        f_dict["disco_num_constraints"] = dist_const["num_constraints"][i]
        f_dict["disco_fraction_observed"] = \
        f_dict["disco_counts"] / f_dict["disco_num_constraints"]

      scores.append(scorer.Predict(f_dict))

    lscores=Table(['chain', 'rindex', 'rnum', 'rname', 'counts', 'packing', 
                   'cb_packing', 'interaction', 'cbeta', 'reduced', 'torsion', 
                   'ss_agreement', 'acc_agreement', 'dist_const', 'QMEAN'],
                   'siisfffffffffff')

    for i, res in enumerate(model.residues):
      lscores.AddRow({'chain' : res.chain.name,
                      'rindex' : res.index,
                      'rnum' : res.number.num,
                      'rname' : res.name,
                      'counts' : data['counts'][i],
                      'packing' : data['packing'][i],
                      'cb_packing' : data['cb_packing'][i],
                      'interaction' : data['interaction'][i],
                      'cbeta' : data['cbeta'][i],
                      'reduced' : data['reduced'][i],
                      'torsion' : data['torsion'][i],
                      'ss_agreement' : data['ss_agreement'][i],
                      'acc_agreement' : data['acc_agreement'][i],
                      'dist_const' :data['dist_const']['disco'][i],
                      'QMEAN' : qmean_scores[i]})

    lscores.comment=LSCORES_TABLE_HEADER

    if assign_bfactors:
      for r,s in zip(model.residues, qmean_scores):
        for a in r.atoms:
          a.b_factor = s

    return LocalResult(model, lscores, local_mqa)



def AssessModelQuality(model, output_dir='.', plots=True, local_scores=True,
                       global_scores=True, table_format='ost', psipred=None, 
                       accpro=None, dc=None, dssp_path=None,
                       assign_bfactors=True, settings=None):
  #hack, that torsion handles get assigned
  processor = conop.HeuristicProcessor(connect=False,peptide_bonds=False,
                                       assign_torsions=True)
  processor.Process(model.handle,False)

  dssp.AssignDSSP(model, extract_burial_status=True, dssp_bin = dssp_path)

  results = []

  if settings == None:
    settings = conf.SwissmodelSettings()

  global_result = None

  if not os.path.exists(output_dir):
    os.makedirs(output_dir)
  plot_dir=os.path.join(output_dir,'plots')
  sliders_dir=os.path.join(plot_dir,'sliders')
  if plots:
    if not os.path.exists(sliders_dir):
      os.makedirs(sliders_dir)
  if global_scores:
    global_result=GlobalResult.Create(model,settings, psipred=psipred, accpro=accpro)
    tab=global_result.score_table
    tab.Save(os.path.join(output_dir,'global_scores.txt'),format=table_format)
    if plots:
      p=global_result.ReferenceSetPlot('QMEAN4')
      p.savefig(os.path.join(plot_dir,'reference_set.png'))
      if psipred!=None and accpro!=None:
        p=global_result.ReferenceSetPlot('QMEAN6')
        p.savefig(os.path.join(plot_dir,'qmean6_reference_set.png'))
      qmean4_scores = [global_result.qmean4, global_result.interaction, global_result.cbeta, global_result.packing, global_result.torsion]
      p=global_result.ZScoreSliders(qmean4_scores)
      p.savefig(os.path.join(plot_dir,'qmean4.png'))
      p.close("all")
      for s in qmean4_scores:
        p=global_result.ZScoreSlider(s)
        p.savefig(os.path.join(sliders_dir,'%s_compact.png'%(sm_conversions[s.name])))
        p=global_result.ZScoreSlider(s, compact=False)
        p.savefig(os.path.join(sliders_dir,'%s.png'%(sm_conversions[s.name])))
        p.close("all")
      if psipred!=None and accpro!=None:
        qmean6_scores = qmean4_scores
        qmean6_scores[0] = global_result.qmean6
        qmean6_scores+=[global_result.ss_agreement,global_result.acc_agreement]
        p=global_result.ZScoreSliders(qmean6_scores)
        p.savefig(os.path.join(plot_dir,'qmean6.png'))
        p.close("all")
        for s in [global_result.qmean6,global_result.ss_agreement,global_result.acc_agreement]:
          p=global_result.ZScoreSlider(s)
          p.savefig(os.path.join(sliders_dir,'%s_compact.png'%(sm_conversions[s.name])))
          p=global_result.ZScoreSlider(s,compact=False)
          p.savefig(os.path.join(sliders_dir,'%s.png'%(sm_conversions[s.name])))
          p.close("all")
    results.append(global_result)
  if local_scores:
    local_result=LocalResult.Create(model,settings,assign_bfactors,psipred=psipred,
                                    accpro=accpro,dc=dc, global_result=global_result)
    tab=local_result.score_table
    tab.Save(os.path.join(output_dir,'local_scores.txt'),format=table_format)
    if plots:
      for ch in model.chains:
        p=local_result.PlotlDDTProfile(chain=ch.name, plot_disco=(dc!=None))
        p.savefig(os.path.join(plot_dir,'local_quality_estimate_%s.png' % (ch.name)))
      p=local_result.PlotlDDTProfile(plot_disco=(dc!=None))
      p.savefig(os.path.join(plot_dir,'local_quality_estimate.png'))
    results.append(local_result)

  return results
