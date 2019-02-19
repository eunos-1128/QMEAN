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
from qmean import mqa_membrane
from qmean import score_calculator
from qmean import reference_set
from qmean import conf
from ost.table import *
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def GenerateEnergyGapPlot(energy, path):
    num_bins = 100
    min_energy = -514128.625
    max_energy = 0.0

    membrane_hist = [
        1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
        0, 0, 0, 0, 0, 0, 0, 2, 0, 1, 4, 0, 1, 0, 0, 1, 0, 2, 1, 0, 1, 0, 0, 1,
        1, 0, 0, 2, 0, 0, 1, 0, 1, 4, 1, 6, 4, 2, 5, 3, 4, 3, 7, 3, 6, 10, 5,
        4, 5, 5, 1, 2, 3, 2, 4, 6, 7, 2, 10, 8, 10, 4, 14, 6, 7, 5, 5, 4, 0, 0,
        0, 0, 0, 0, 0, 1
    ]

    soluble_hist = [
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 2, 15, 179
    ]

    bin_edges = [
        -514128.625, -508987.33875, -503846.0525, -498704.76625, -493563.48,
        -488422.19375, -483280.9075, -478139.62125, -472998.335, -467857.04875,
        -462715.7625, -457574.47625, -452433.19, -447291.90375, -442150.6175,
        -437009.33125, -431868.045, -426726.75875, -421585.4725, -416444.18625,
        -411302.9, -406161.61375, -401020.3275, -395879.04125, -390737.755,
        -385596.46875, -380455.1825, -375313.89625, -370172.61, -365031.32375,
        -359890.0375, -354748.75125, -349607.465, -344466.17875, -339324.8925,
        -334183.60625, -329042.32, -323901.03375, -318759.7475, -313618.46125,
        -308477.175, -303335.88875, -298194.6025, -293053.31625, -287912.03,
        -282770.74375, -277629.4575, -272488.17125, -267346.885, -262205.59875,
        -257064.3125, -251923.02625, -246781.74, -241640.45375, -236499.1675,
        -231357.88125, -226216.595, -221075.30875, -215934.0225, -210792.73625,
        -205651.45, -200510.16375, -195368.8775, -190227.59125, -185086.305,
        -179945.01875, -174803.7325, -169662.44625, -164521.16, -159379.87375,
        -154238.5875, -149097.30125, -143956.015, -138814.72875, -133673.4425,
        -128532.15625, -123390.87, -118249.58375, -113108.2975, -107967.01125,
        -102825.725, -97684.43875, -92543.1525, -87401.86625, -82260.58,
        -77119.29375, -71978.0075, -66836.72125, -61695.435, -56554.14875,
        -51412.8625, -46271.57625, -41130.29, -35989.00375, -30847.7175,
        -25706.43125, -20565.145, -15423.85875, -10282.5725, -5141.28625
    ]

    bar_width = (max_energy - min_energy) / num_bins

    fig, ax = plt.subplots()

    ax.bar(
        bin_edges,
        membrane_hist,
        width=bar_width,
        color=(204.0 / 255, 0, 0),
        label='Membrane Structures',
        zorder=1)
    ax.bar(
        bin_edges,
        soluble_hist,
        width=bar_width,
        color=(0.0, 147.0 / 255, 217.0 / 255),
        label='Soluble Structures',
        zorder=1)

    if energy < -250000:
        if energy < -510000:
            text = "Pseudo energy of %.f has been\n capped at %.f" % (
                energy, -510000.00)
            # for the usage later on
            energy = -510000
            ax.text(
                -480000,
                80,
                text,
                color='black',
                bbox=dict(
                    facecolor="#ffb2b2",
                    alpha=1.0,
                    edgecolor='black',
                    boxstyle='round,pad=1',
                    fill=True),
                zorder=3)

        ax.plot([energy, energy], [0.0, 150], 'black', linewidth=4.0, zorder=2)
    else:
        ax.plot(
            [energy, energy], [0.0, 177.9], 'black', linewidth=4.0, zorder=2)

    # Annotation for the model specific energy
    text_coord = None
    c_style = None
    if energy < -250000:
        text_coord = (energy + 10000, 120)
        c_style = "arc3,rad=-0.3"
    else:
        text_coord = (energy - 100000, 120)
        c_style = "arc3,rad=0.3"

    ax.annotate(
        'Your Model',
        xy=(energy, 100),
        xytext=text_coord,
        fontsize='large',
        arrowprops=dict(
            arrowstyle="simple", fc="0.6", ec="none", connectionstyle=c_style))

    # info box
    if (energy < -30000):
        text = "Pseudo energy output of membrane\n"
        text += "locating algorithm. Your model\n"
        text += "is within the expected range of\n"
        text += "a transmembrane structure.\n"

        ax.text(
            -480000,
            20,
            text,
            color='black',
            bbox=dict(
                facecolor="#b2d8b2",
                alpha=1.0,
                edgecolor='black',
                boxstyle='round,pad=1'),
            zorder=3)

    else:
        text = "Pseudo energy output of membrane\n"
        text += "locating algorithm. Are you sure\n"
        text += "that your model represents a\n"
        text += "transmembrane structure?\n"
        text += "It might be of very low quality or\n"
        text += "of wrong quarternary state making\n"
        text += "it difficult to locate the membrane!"

        ax.text(
            -480000,
            20,
            text,
            color='black',
            bbox=dict(
                facecolor="#ffb2b2",
                alpha=1.0,
                edgecolor='black',
                boxstyle='round,pad=1',
                fill=True),
            zorder=3)

    # set limits, produce legends, axis description...
    plt.xlim((min_energy, max_energy))
    plt.legend(loc=2, frameon=False, fontsize='large')
    plt.ylabel('Structures', fontsize='large')
    plt.xlabel('Membrane Insertion Energy', fontsize='large')

    xtick_positions = [-500000, -400000, -300000, -200000, -100000, 0]
    xtick_labels = ['-500k', '-400k', '-300k', '-200k', '-100k', '0']

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.set_xticks(xtick_positions, xtick_labels)

    fig.savefig(path)





LSCORES_TABLE_HEADER='''\
This file contains local scores calculated by the QMEAN scoring function.

References:

  Studer, G., Biasini, M. and Schwede, T. (2014). "Assessing the local 
  structural quality of transmembrane protein models using statistical 
  potentials (QMEANBrane)" Bioinformatics, 30(17): i505-i511.

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

  QMEAN         Predicted local quality. A value describing the expected local
                similarity to the target structure with a range of [0,1]. 
'''

class LocalMembraneResult:
  def __init__(self,model,data):
    self.model = model
    self.score_table = data

  def PlotlDDTProfile(self, chain=None):

    from matplotlib import pyplot
    pyplot.clf()
    pyplot.figure(figsize=[8,6],dpi=80)
    pyplot.ylim((-0.1,1.1))
    pyplot.xlabel('Residue Number',size='large')
    pyplot.ylabel('Predicted Local Similarity to Target',size='large')

#    color_scheme = ['#D55E00','#0072B2','#F0E442','#009E73','#56B4E9','#E69F00','#999999','#000000']
    color_scheme = ['#ff420e','#004586','#ffd320','#578d1c','#7e0021']

    chain_idx = self.score_table.GetColIndex('chain')
    rnum_idx = self.score_table.GetColIndex('rnum')
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
  def Create(model, settings, assign_bfactors, membrane_query = None, interface_query = None, mem_param = None, psipred=None, accpro=None):

    pot_membrane = PotentialContainer.Load(settings.local_potentials_membrane)
    pot_soluble = PotentialContainer.Load(settings.local_potentials_soluble)
    scorer_membrane = score_calculator.LocalScorer.Load(settings.local_scorer_membrane)
    scorer_soluble = score_calculator.LocalScorer.Load(settings.local_scorer_soluble)
    local_mqa = mqa_membrane.MembraneScores(model, model, pot_soluble, pot_membrane, smooth_std=5.0, psipred=psipred, 
                                            accpro=accpro, membrane_query = membrane_query, 
                                            interface_query = interface_query,mem_param = mem_param)


    features = ['interaction','cbeta','packing','torsion','exposed']
    if psipred != None:
      features.append('ss_agreement')
    if accpro != None:
      features.append('acc_agreement')

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
        scores.append(min(max(scorer_soluble.GetLocalScore(ss, r.one_letter_code, residue_data),0.0),1.0))

      else:
        residue_data=dict()
        for f in features:
          if f in ['ss_agreement','acc_agreement']:
            continue
          residue_data[f] = data[f][i]
        if membrane_states[i] == 1:
          scores.append(min(max(scorer_membrane.GetLocalScore('membrane', r.one_letter_code, residue_data),0.0),1.0))
        elif membrane_states[i] == 2:
          scores.append(min(max(scorer_membrane.GetLocalScore('interface', r.one_letter_code, residue_data),0.0),1.0))
        else:
          raise RuntimeError("Invalid membrane state encountered")
          

    data['QMEAN'] = scores


    if assign_bfactors:
      for r,s in zip(model.residues,data['QMEAN']):
        for a in r.atoms:
          a.b_factor = s


    #check whether ss_agreement and acc_agreement are set, otherwise we just
    #add NaN values
    if 'ss_agreement' not in features:
      data['ss_agreement'] = [float('NaN')] * len(data[features[0]])
    if 'acc_agreement' not in features:
      data['acc_agreement'] = [float('NaN')] * len(data[features[0]])



    lscores=Table(['chain', 'rindex', 'rnum', 'rname', 'membrane','all_atom',
                   'cbeta', 'solvation', 'torsion', 'exposed',
                   'ss_agreement', 'acc_agreement', 'QMEAN'],
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
                      'QMEAN' : data['QMEAN'][i]})

    lscores.comment=LSCORES_TABLE_HEADER
    return LocalMembraneResult(model,lscores)


def AssessMembraneModelQuality(model, mem_param = None, output_dir='.', 
                               plots=True, table_format='ost', psipred=None, 
                               accpro=None, assign_bfactors=True, settings=None):

  if settings == None:
    settings = conf.MembraneSettings()

  membrane_query = None
  interface_query = None

  #hack, that torsion handles get assigned
  processor = conop.HeuristicProcessor(connect=False,peptide_bonds=False,assign_torsions=True)
  processor.Process(model.handle,False)

  if not os.path.exists(output_dir):
    os.makedirs(output_dir)
  plot_dir=os.path.join(output_dir,'plots')
  if plots:
    if not os.path.exists(plot_dir):
      os.makedirs(plot_dir)
  local_result=LocalMembraneResult.Create(model,settings,assign_bfactors,
                                          membrane_query, interface_query,
                                          mem_param, psipred=psipred, 
                                          accpro=accpro)
  tab=local_result.score_table
  tab.Save(os.path.join(output_dir,'local_scores.txt'),format=table_format)
  if plots:
    for ch in model.chains:
      p=local_result.PlotlDDTProfile(chain=ch.name)
      p.savefig(os.path.join(plot_dir,'local_quality_estimate_%s.png' % (ch.name)))
    p=local_result.PlotlDDTProfile()
    p.savefig(os.path.join(plot_dir,'local_quality_estimate.png'))

  return local_result

