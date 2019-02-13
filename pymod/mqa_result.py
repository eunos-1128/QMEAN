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

from qmean import mqa
from qmean import score_calculator
from qmean import reference_set
from qmean import conf
from qmean import PotentialContainer
from ost.table import Table
from ost import mol
import matplotlib.pyplot as plt

class QMEANScorer(object):

  def __init__(self, model, psipred = None, accpro = None, dc = None, 
               settings = conf.SwissmodelSettings(), use_nn = True):

    self._model = model.Select("peptide=true")
    # We rely on residue numbers to uniquely identify residues in a chain.
    # There's no guarantee that they're are unique, so we check manually
    for c in self._model.chains:
      residue_numbers = set()
      for r in c.residues:
        num = r.GetNumber()
        # no insertion codes allowed
        if len(num.GetInsCode().strip('\0')) > 0:
          err = "QMEAN cannot deal with insertion codes in residue numbers. "
          err += "Problematic residue: " + str(r)
          raise RuntimeError(err)

        residue_numbers.add(r.GetNumber().GetNum())

      if len(residue_numbers) != len(c.residues):
        err = "The residue numbers within one chain must be unique. "
        err += "This is not the case for chain " + str(c)
        raise RuntimeError(err)

    self._psipred = psipred
    self._accpro = accpro
    self._dc = dc
    self._settings = settings
    self.use_nn = use_nn
    self._local_mqa = None
    self._global_mqa = None
    self._linear_global_scorer = None
    self._linear_local_scorer = None
    self._classic_reference_set = None
    self._nn_scorer = None
    self._qmean4_score = None
    self._qmean4_z_score = None
    self._qmean6_score = None
    self._qmean6_z_score = None
    self._local_scores = None
    self._avg_local_score = None


  @property
  def local_mqa(self):
    if self._local_mqa == None:
      pot = PotentialContainer.Load(self._settings.local_potentials)
      self._local_mqa = mqa.Scores(self._model, self._model, pot, 
                                   smooth_std=5.0, psipred=self._psipred, 
                                   accpro=self._accpro, dc=self._dc)
    return self._local_mqa


  @property
  def global_mqa(self):
    if self._global_mqa == None:
      pot = PotentialContainer.Load(self._settings.global_potentials)
      self._global_mqa = mqa.Scores(self._model, self._model, pot, 
                                    psipred=self._psipred, accpro=self._accpro)
    return self._global_mqa


  @property
  def linear_global_scorer(self):
    if self._linear_global_scorer == None:
      self._linear_global_scorer = \
      score_calculator.GlobalScorer.Load(self._settings.global_scorer)
    return self._linear_global_scorer


  @property
  def linear_local_scorer(self):
    if self._linear_local_scorer == None:
      self._linear_local_scorer = \
      score_calculator.LocalScorer.Load(self._settings.linear_local_scorer)
    return self._linear_local_scorer


  @property
  def classic_reference_set(self):
    if self._classic_reference_set == None:
      t = Table.Load(self._settings.reference_tab)
      self._classic_reference_set = reference_set.ReferenceSet(t)
    return self._classic_reference_set


  @property
  def nn_scorer(self):
    if self._nn_scorer == None:
      self._nn_scorer = score_calculator.NNScorer(self._settings.local_scorer)
    return self._nn_scorer
  

  @property
  def qmean4_score(self):
    if self._qmean4_score == None:
      data = self.global_mqa.GetAVGData(['interaction','cbeta',
                                         'packing','torsion'])
      s = self.linear_global_scorer.GetGlobalScore('soluble', data)
      self._qmean4_score = min(max(s, 0.0), 1.0)
    return self._qmean4_score


  @property
  def qmean4_z_score(self):
    if self._qmean4_z_score == None:
      ref_set = self.classic_reference_set
      s = ref_set.ZScoreFromNormScore('QMEAN4', self._model.residue_count,
                                      self.qmean4_score)
      self._qmean4_z_score = s
    return self._qmean4_z_score


  @property
  def qmean6_score(self):
    if self._qmean6_score == None:
      if self._psipred == None or self._accpro == None:
        raise RuntimeError("Require PSIPRED and ACCPRO information for QMEAN6")
      data = self.global_mqa.GetAVGData(['interaction','cbeta',
                                         'packing','torsion',
                                         'ss_agreement', 'acc_agreement'])
      s = self.linear_global_scorer.GetGlobalScore('soluble', data)
      self._qmean6_score = min(max(s, 0.0), 1.0)
    return self._qmean6_score


  @property
  def qmean6_z_score(self):
    if self._qmean6_z_score == None:
      ref_set = self.classic_reference_set
      s = ref_set.ZScoreFromNormScore('QMEAN6', self._model.residue_count,
                                      self.qmean6_score)
      self._qmean6_z_score = s
    return self._qmean6_z_score


  @property
  def local_scores(self):
    if self._local_scores == None:
      self._local_scores = dict()

      if self.use_nn:
        local_features = ["counts", "packing", "cb_packing", "exposed", 
                          "torsion", "reduced", "interaction", "cbeta", 
                          "ss_agreement", "acc_agreement", "dist_const", 
                          "clash"]
        avg_features = ["torsion", "reduced", "interaction", "cbeta", 
                        "packing", "cb_packing", "ss_agreement", 
                        "acc_agreement"]
        local_data = self.local_mqa.GetLocalData(local_features)
        avg_data = self.local_mqa.GetAVGData(avg_features)
        scorer = self.nn_scorer
      
        for i, r in enumerate(self._model.residues):
          f_dict = dict()

          # do the avg features
          for f in avg_features:
            f_dict["avg_" + f] = avg_data[f]

          # do the local features
          for f in local_features:
            # dist const requires special treatment
            if f == "dist_const":
              continue
            f_dict[f] = local_data[f][i]

          # do dist_const
          dc = local_data["dist_const"]
          s = dc["disco"][i]
          if s == s:
            f_dict["disco"] = s
            f_dict["disco_counts"] = dc["counts"][i]
            f_dict["disco_avg_num_clusters"] = dc["avg_num_clusters"][i]
            f_dict["disco_avg_max_seqsim"] = dc["avg_max_seqsim"][i]
            f_dict["disco_avg_max_seqid"] = dc["avg_max_seqid"][i]
            f_dict["disco_avg_variance"] = dc["avg_variance"][i]
            f_dict["disco_num_constraints"] = dc["num_constraints"][i]
            f_dict["disco_fraction_observed"] = dc["fraction_observed"][i]

          # predict score and make sure that the range is in [0.0, 1.0]
          s = scorer.GetScore(f_dict, r.one_letter_code)
          s = min(1.0, max(0.0, s))

          chain_name = r.GetChain().GetName()
          if chain_name not in self._local_scores:
            self._local_scores[chain_name] = dict()
          self._local_scores[chain_name][r.GetNumber().GetNum()] = s
      else:
        ##############
        # DEPRECATED #
        ##############
        features = ["interaction", "cbeta", "packing", "torsion", "exposed",
                    "ss_agreement", "acc_agreement"]
        data = self.local_mqa.GetLocalData(features)
        scorer = self.linear_local_scorer
        dssp_ss = self.local_mqa.GetDSSPSS()
  
        for i ,r in enumerate(self._model.residues):
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
          s = scorer.GetLocalScore(ss, r.one_letter_code, residue_data)
          s = min(1.0, max(0.0, s))
  
          chain_name = r.GetChain().GetName()
          if chain_name not in self._local_scores:
            self._local_scores[chain_name] = dict()
          self._local_scores[chain_name][r.GetNumber().GetNum()] = s

    return self._local_scores


  @property
  def avg_local_score(self):
    if self._avg_local_score == None:
      self._avg_local_score = 0.0
      for chain_name, residue_scores in self._local_scores.iteritems():
        for rnum, s in residue_scores.iteritems():
          self._avg_local_score += s
      self._avg_local_score /= self._model.residue_count
    return self._avg_local_score


  def QMEAN4SliderPlot(self, out_path):
    data = self.global_mqa.GetAVGData(['interaction','cbeta',
                                       'packing','torsion'])
    norm_scores = [self.qmean4_score, data['interaction'], data['cbeta'],
                   data['packing'], data['torsion']]
    # naming in the ref set
    names = ['QMEAN4', 'interaction', 'cbeta', 'packing', 'torsion']
    # naming in the plot
    labels = ['QMEAN4', 'All Atom', 'CBeta', 'Solvation', 'Torsion']
    self.classic_reference_set.DumpZScoreSliders(out_path, norm_scores, names, 
                                                 self._model.residue_count, 
                                                 score_labels = labels)


  def QMEAN6SliderPlot(self, out_path):
    data = self.global_mqa.GetAVGData(['interaction','cbeta',
                                       'packing','torsion',
                                       'ss_agreement', 'acc_agreement'])
    norm_scores = [self.qmean6_score, data['interaction'], data['cbeta'],
                   data['packing'], data['torsion'], data['ss_agreement'],
                   data['acc_agreement']]
    # naming in the ref set
    names = ['QMEAN6', 'interaction', 'cbeta', 'packing', 'torsion', 
             'ss_agreement', 'acc_agreement']
    # naming in the plot
    labels = ['QMEAN6', 'All Atom', 'CBeta', 'Solvation', 'Torsion',
              'SS Agree', 'Acc Agree']
    self.classic_reference_set.DumpZScoreSliders(out_path, norm_scores, names, 
                                                 self._model.residue_count, 
                                                 score_labels = labels)


  def QMEAN4ReferencePlot(self, out_path):
    self.classic_reference_set.DumpReferenceSetPlot(out_path, 
                                                    self._model.residue_count, 
                                                    self.qmean4_score, 'QMEAN4')


  def QMEAN6ReferencePlot(self, out_path):
    self.classic_reference_set.DumpReferenceSetPlot(out_path, 
                                                    self._model.residue_count, 
                                                    self.qmean6_score, 'QMEAN6')


  def LocalProfilePlot(self, out_path, chain = None,
                       colors = ['#ff420e', '#004586', '#ffd320', 
                                 '#578d1c', '#7e0021']):
    if chain != None and chain not in self.local_scores:
      err = "Can not plot local score profile for chain " + chain
      err += ". Did not calculate any local scores for it."
      raise RuntimeError(err)

    plt.clf()
    plt.figure(figsize=[8,6],dpi=80)
    plt.ylim((-0.1,1.1))
    plt.xlabel('Residue Number',size='large')
    plt.ylabel('Predicted Local Similarity to Target', size='large')

    if chain == None:
      chains = self.local_scores.keys()
      chains.sort()
      plt.title('Local Quality Estimate',size='x-large')
      for ch_idx, ch in enumerate(chains):
        data = list()
        for rnum, score in self.local_scores[ch].iteritems():
          data.append((rnum, score))
        data.sort()
        res_num = [x[0] for x in data]
        scores = [x[1] for x in data]
        plt.plot(res_num, scores, color = colors[ch_idx%len(colors)], 
                 linewidth = 2.0)         
    else:
      plt.title('Local Quality Estimate: Chain %s'%(chain), size = 'x-large')
      data = list()
      for rnum, score in self.local_scores[chain].iteritems():
        data.append((rnum, score))
      data.sort()
      res_num = [x[0] for x in data]
      scores = [x[1] for x in data]
      plt.plot(res_num, scores, color = colors[0], linewidth = 2.0) 

    plt.savefig(out_path)

  def GetModel(self):
    return self._model

  def AssignModelBFactors(self):
    for chain_name, data in self.local_scores.iteritems():
      chain = self._model.FindChain(chain_name)
      for num, score in data.iteritems():
        res = chain.FindResidue(mol.ResNum(num))
        for a in res.atoms:
          a.b_factor = score
        
    
