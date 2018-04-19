import bisect
import math
from ost import stutil


class ReferenceSet:
  """
  Holds the reference QMEAN values calculated for a large set of X-ray
  structures.
  """
  def __init__(self, tab):
    self.data=tab
    # we will have to perform lookups based on protein size. Doing this with a
    # binary search really helps. So let's sort by size
    self.data.Sort('num_residues', '-')

    self.nice_titles={'interaction': 'Normalised All Atom Pairwise Score',
            'cbeta': 'Normalized CBeta Pairwise Score',
            'torsion': 'Normalized Torsion Score',
            'packing': 'Normalized Solvation Score',
            'ss_agreement': 'Normalized SSE Agreement Score',
            'acc_agreement': 'Normalized ACC Agreement Score',
            'reduced':'Normalized Reduced Score',
            'QMEAN4': 'Normalized QMEAN4 Score',
            'QMEAN6': 'Normalized QMEAN6 Score'}

  def ZScoreFromNormScore(self, score_name, num_residues, norm_score):
    """
    Calculates the ZScore from the normalized QMEAN score.

    :param score_name: Name of the term. Must be one of "all atom", "solvation",
       "cbeta", "torsion", "sse", "acc". More will likely come in the future.
    :param norm_score: The normalized score
    :param num_residues: Number of residues of the model structure.
    :raises: :exc:`ValueError` if *score_name* is not one of the accepted
       values.
    """

    inverted_scores = set(['QMEAN4', 'QMEAN6', 'ss_agreement', 
                           'acc_agreement'])
    if score_name not in self.nice_titles:
      raise ValueError('unknown score term "%s"' % score_name)
    num_residues=max(40, min(num_residues, 600))
    first=bisect.bisect_left(self.data['num_residues'], int(num_residues*0.9))
    last=bisect.bisect_right(self.data['num_residues'], int(num_residues*1.1))
    score_values=list(self.data[score_name])[first:last]
    mean=stutil.Mean(score_values)
    
    sign = 1
    if score_name in inverted_scores:
      sign = -1
    
    # we can't use stutil.StdDev directly, since we have to use the n-1
    # convention here... That's what R is doing in the "sd" function.
    diff=0

    for val in score_values:
      diff+=(val-mean)**2
    stdev=math.sqrt(diff/(len(score_values)-1))
    # I know, usually Zscores are written as (score-mean)/stdev. QMean uses the
    # convention that high values are good while small values are bad, so we
    # have to invert the sign...
    return sign*(mean-norm_score)/stdev

  def DumpReferenceSetPlot(self, out_path, num_residues, qmean_norm_score, 
                           score_name):
    """
    Plots and dumps the qmean_norm_score in relation to QMEAN scores for a 
    non-redundant set of X-ray structures.
    """
    if score_name not in self.nice_titles:
      raise ValueError('invalid score: '+score_name)
    col_index=self.data.GetColIndex('ZScore_'+score_name)
    cutoff=1000
    def _SelectValues(table, row):
      return cutoff>abs(row[col_index])
    style='.'
    self.data.Plot('num_residues', score_name, style=style, color=(0.8,0.8, 0.8),
                   legend='|Z-score|>2',
                   plot_if=_SelectValues)
    cutoff=2
    self.data.Plot('num_residues', score_name, style=style,
                   color=(0.6,0.6, 0.6),
                   legend='1<|Z-score|<2',
                   plot_if=_SelectValues, clear=False)
    cutoff=1
    p=self.data.Plot('num_residues', score_name, x_range=(0, 600), style=style,
                     clear=False, legend='|Z-score|<1',
                     y_range=(-0.5, 1.3), y_title=self.nice_titles[score_name],
                     x_title='Protein Size (Residues)', color=(0.3,0.3,0.3),
                     title='Comparison with Non-redundant Set of PDB Structures',
                     plot_if=_SelectValues)
    p.plot([min(num_residues,595)], [qmean_norm_score], 'r*', label='model',
           markeredgecolor='white', markersize=10)
    p.legend(loc=4, markerscale=2,numpoints=1, frameon=False)

    #these things have already been set, in the table plot class
    #Nevertheless we overwrite it to make it consistent with the 
    #style of other plots
    p.title('Comparison with Non-redundant Set of PDB Structures',
             size='x-large', fontweight = 'normal')
    p.xlabel('Protein Size (Residues)', size='large',fontweight='normal')
    p.ylabel(self.nice_titles[score_name],size='large',fontweight='normal')

    if(num_residues>600):
      text = 'Number of residues (%i) has been\n'%(num_residues)
      text += 'capped at 600 for visualisation\nand Z-score calculation!'
      p.text(35,-0.365,text,color='black',bbox=dict(facecolor='red', alpha=0.3,
                                                    edgecolor='black',
                                                    boxstyle='round,pad=1'))

    p.savefig(out_path)

  def DumpZScoreSliders(self, out_path, qmean_norm_scores, score_names, 
                        num_residues, score_labels=None, **kwargs):
    """
    Plots and dumps QMEAN scores as sliders.
    """
    import matplotlib as mpl
    from matplotlib import pyplot
    pyplot.clf()

    if len(qmean_norm_scores)!=len(score_names):
      raise ValueError('Number of scores and score names is not consistent!')

    if score_labels!=None:
      if len(qmean_norm_scores)!=len(score_names):
        raise ValueError('Number of scores and score labels is not consistent!')
    else:
      score_labels = score_names
    
    compact=kwargs.get('compact', False)
    if compact:
      if len(qmean_norm_scores)>1:
        raise ValueError('Only one term is allowed in compact mode!')


    if not compact:
      full_height=1+0.3*len(score_names)
      fig = pyplot.figure(figsize=[6, full_height])
      pyplot.suptitle('Z-score', size='x-large')
      slider_height=(1.0/(2+len(score_names)))*0.8
      slider_space=1.0/(2+len(score_names))
      #slider_height=0.2
      top_offset=1.0/(2+len(score_names))
      slider_start=0.2
      slider_width=0.65
    else:
      fig = pyplot.figure(figsize=[3,.3])
      slider_height=1.0
      top_offset=0.00
      slider_start=0.00
      slider_width=1.00

    cdict = {'red':  [(0.0,1.0,1.0),
                      (0.75,1.0,1.0),
                      (1.0,0.3,0.3)],
             'green':[(0.0,0.3,0.3),
                      (0.75,1.0,1.0),
                      (1.0,0.3,0.3)],
             'blue': [(0.0,0.3,0.3),
                      (0.75,1.0,1.0),
                      (1.0,1.0,1.0)]}

    bwr_gradient = mpl.colors.LinearSegmentedColormap('bwr',cdict)

    slider=1
    for score_name, score, score_label in zip(score_names, qmean_norm_scores, score_labels):
      z_score=self.ZScoreFromNormScore(score_name,num_residues,score)
      if compact:
        y_pos = 0.0
      else:
        y_pos=(1-top_offset-slider*slider_space)
        fig.text(slider_start-0.02, y_pos+(slider_space-slider_height), score_label,
               fontsize=14, horizontalalignment='right',
               verticalalignment='baseline')
        fig.text(slider_start+slider_width+0.1, y_pos+(slider_space-slider_height),
               '%.2f' % z_score, fontsize=14, horizontalalignment='right',
               verticalalignment='baseline')
      
      ax= fig.add_axes([slider_start, y_pos, slider_width, slider_height],
                       xlim=[-6, 2], ylim=[0, 1], autoscalex_on=False,
                       autoscaley_on=False)
      
      
      bounded_zscore=max(-5.94, min(z_score, 2.0))
      transf_zscore=(bounded_zscore+6)/8.0

      ax.plot([transf_zscore, transf_zscore], [0, 1],
              linewidth=3, color='black')

      ax.set_frame_on(False)
      norm = mpl.colors.Normalize(vmin=-6, vmax=2)
      cb = mpl.colorbar.ColorbarBase(ax, cmap=bwr_gradient,
                                    norm=norm,
                                    orientation='horizontal',
                                    ticks=[-6, -5, -4, -3, -2, -1, 0, 1, 2])
      if slider<len(score_names):
        ax.set_xticklabels([])

      if compact:
        ax.set_xticklabels([])
      slider+=1

    pyplot.savefig(out_path)
