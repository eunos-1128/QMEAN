import os
from scipy.ndimage import gaussian_filter
from ost.table import Table
import matplotlib.pyplot as plt
import numpy as np


#Load the data produced by QMEANBrane
tab_no_shift = Table.Load(os.path.join('original_hhblits_alignment',
                                       'local_scores.txt'))
tab_correct_shift = Table.Load(os.path.join('shift_in_front_helix_four',
                                            'local_scores.txt'))
tab_wrong_shift = Table.Load(os.path.join('shift_towards_cter',
                                          'local_scores.txt'))
tab_shift_into_middle = Table.Load(os.path.join('shift_into_middle',
                                                'local_scores.txt'))

score_idx = tab_no_shift.GetColIndex('QMEAN')
rnum_idx = tab_no_shift.GetColIndex('rnum')

rnums_no_shift = list()
scores_no_shift = list()
for r in tab_no_shift.rows:
  rnums_no_shift.append(r[rnum_idx])
  scores_no_shift.append(r[score_idx])

scores_correct_shift = dict()
for r in tab_correct_shift.rows:
  scores_correct_shift[r[rnum_idx]] = r[score_idx]

scores_wrong_shift = dict()
for r in tab_wrong_shift.rows:
  scores_wrong_shift[r[rnum_idx]] = r[score_idx]

scores_shift_into_middle = dict()
for r in tab_shift_into_middle.rows:
  scores_shift_into_middle[r[rnum_idx]] = r[score_idx]

rnums = list()

#calculate the differences of the scores originating from models built
#with modified alignments relative to the model built with the original 
#hhblits alignment.

differences_correct_shift = list()
differences_wrong_shift = list()
differences_shift_into_middle = list()

for n,s in zip(rnums_no_shift,scores_no_shift):
  if n not in scores_correct_shift:
    continue
  if n not in scores_wrong_shift:
    continue
  if n not in scores_shift_into_middle:
    continue
  rnums.append(n)
  differences_correct_shift.append(scores_correct_shift[n]-s)
  differences_wrong_shift.append(scores_wrong_shift[n]-s)
  differences_shift_into_middle.append(scores_shift_into_middle[n]-s)

#do some slight smoothing
differences_correct_shift = gaussian_filter(differences_correct_shift,1)
differences_wrong_shift = gaussian_filter(differences_wrong_shift,1)
differences_shift_into_middle = gaussian_filter(differences_shift_into_middle,1)

#membrane spanners define the transmembrane helices
membrane_spanners = list()
membrane_spanners.append([11,37])
membrane_spanners.append([43,74])
membrane_spanners.append([86,121])
membrane_spanners.append([175,193])
membrane_spanners.append([198,222])
membrane_spanners.append([246,271])
membrane_spanners.append([279,307])
membrane_spanners.append([341,374])
membrane_spanners.append([382,402])
membrane_spanners.append([408,434])
membrane_spanners.append([456,479])
membrane_spanners.append([497,526])

#following lists contain information regarding insertions and deletions
#in the four different alignments

#initial alignment
insertions = list()
insertions.append((1,5))
insertions.append((127,127))
insertions.append((132,134))
insertions.append((140,144))
insertions.append((162,163))
insertions.append((204,206))
insertions.append((308,311))
insertions.append((378,380))
insertions.append((406,407))
insertions.append((487,492))
insertions.append((529,543))

deletions = list()
deletions.append((85,86))
deletions.append((166,167))
deletions.append((223,224))
deletions.append((294,295))
deletions.append((305,306))
deletions.append((319,320))
deletions.append((321,322))

matches = list()
matches.append((5,84))
matches.append((87,126))
matches.append((128,131))
matches.append((135,139))
matches.append((145,161))
matches.append((164,165))
matches.append((168,203))
matches.append((207,222))
matches.append((225,293))
matches.append((296,304))
matches.append((307,307))
matches.append((312,318))
matches.append((323,377))
matches.append((381,405))
matches.append((408,486))
matches.append((493,528))

#correct alignment
insertions_correct = list()
insertions_correct.append((1,5))
insertions_correct.append((127,127))
insertions_correct.append((132,134))
insertions_correct.append((140,144))
insertions_correct.append((162,163))
insertions_correct.append((171,174))
insertions_correct.append((308,311))
insertions_correct.append((378,380))
insertions_correct.append((406,407))
insertions_correct.append((487,492))
insertions_correct.append((529,543))

deletions_correct = list()
deletions_correct.append((85,86))
deletions_correct.append((166,167))
deletions_correct.append((197,198))
deletions_correct.append((223,224))
deletions_correct.append((294,295))
deletions_correct.append((305,306))
deletions_correct.append((319,320))
deletions_correct.append((321,322))

matches_correct = list()
matches_correct.append((5,84))
matches_correct.append((87,126))
matches_correct.append((128,131))
matches_correct.append((135,139))
matches_correct.append((145,161))
matches_correct.append((164,165))
matches_correct.append((168,170))
matches_correct.append((175,196))
matches_correct.append((199,222))
matches_correct.append((225,293))
matches_correct.append((296,304))
matches_correct.append((307,307))
matches_correct.append((312,318))
matches_correct.append((323,377))
matches_correct.append((381,405))
matches_correct.append((408,486))
matches_correct.append((493,528))

#shifted into middle
insertions_middle = list()
insertions_middle.append((1,5))
insertions_middle.append((127,127))
insertions_middle.append((132,134))
insertions_middle.append((140,144))
insertions_middle.append((162,163))
insertions_middle.append((190,190))
insertions_middle.append((191,192))
insertions_middle.append((194,195))
insertions_middle.append((308,311))
insertions_middle.append((378,380))
insertions_middle.append((406,407))
insertions_middle.append((487,492))
insertions_middle.append((529,543))

deletions_middle = list()
deletions_middle.append((85,86))
deletions_middle.append((166,167))
deletions_middle.append((189,190))
deletions_middle.append((190,191))
deletions_middle.append((223,224))
deletions_middle.append((294,295))
deletions_middle.append((305,306))
deletions_middle.append((319,320))
deletions_middle.append((321,322))

matches_middle = list()
matches_middle.append((5,84))
matches_middle.append((87,126))
matches_middle.append((128,131))
matches_middle.append((135,139))
matches_middle.append((145,161))
matches_middle.append((164,165))
matches_middle.append((168,188))
matches_middle.append((193,193))
matches_middle.append((196,222))
matches_middle.append((225,293))
matches_middle.append((296,304))
matches_middle.append((307,307))
matches_middle.append((312,318))
matches_middle.append((323,377))
matches_middle.append((381,405))
matches_middle.append((408,486))
matches_middle.append((493,528))

#shifted towards C-Terminus
insertions_cter = list()
insertions_cter.append((1,5))
insertions_cter.append((127,127))
insertions_cter.append((132,134))
insertions_cter.append((140,144))
insertions_cter.append((162,163))
insertions_cter.append((308,311))
insertions_cter.append((378,380))
insertions_cter.append((406,407))
insertions_cter.append((487,492))
insertions_cter.append((529,543))

deletions_cter = list()
deletions_cter.append((85,86))
deletions_cter.append((166,167))
deletions_cter.append((223,224))
deletions_cter.append((294,295))
deletions_cter.append((305,306))
deletions_cter.append((319,320))
deletions_cter.append((321,322))

matches_cter = list()
matches_cter.append((5,84))
matches_cter.append((87,126))
matches_cter.append((128,131))
matches_cter.append((135,139))
matches_cter.append((145,161))
matches_cter.append((164,165))
matches_cter.append((168,222))
matches_cter.append((225,293))
matches_cter.append((296,304))
matches_cter.append((307,307))
matches_cter.append((312,318))
matches_cter.append((323,377))
matches_cter.append((381,405))
matches_cter.append((408,486))
matches_cter.append((493,528))

#define some colors and plot the differences in scores
red = (204.0/255,0.0/255,0.0/255)
blue = (0.0/255,147.0/255,217.0/255)
orange = (255.0/255,141.0/255,0.0)

plt.plot(rnums,differences_shift_into_middle,color=blue)
plt.plot(rnums,differences_wrong_shift,color=orange)
plt.plot(rnums,differences_correct_shift,color=red)

#do spannners for initial alignment
for ins in insertions:
  plt.axvspan(ins[0],ins[1],ymin=0.975,ymax=1.0,
              color='k',edgecolor='none')
for d in deletions:
  plt.axvspan(d[0],d[1],ymin=0.975,ymax=1.0,
              color='white',edgecolor='none')
for m in matches:
  plt.axvspan(m[0],m[1],ymin=0.975,ymax=1.0,
              color='k',alpha=0.4,edgecolor='none')

#do spanners for correct alignment
for ins in insertions_correct:
  plt.axvspan(ins[0],ins[1],ymin=0.95,ymax=0.975,
              color='k',edgecolor='none')
for d in deletions_correct:
  plt.axvspan(d[0],d[1],ymin=0.95,ymax=0.975,
              color='white',edgecolor='none')
for m in matches_correct:
  plt.axvspan(m[0],m[1],ymin=0.95,ymax=0.975,
              color=red,alpha=0.4,edgecolor='none')

#do spanners for shift into middle
for ins in insertions_middle:
  plt.axvspan(ins[0],ins[1],ymin=0.925,ymax=0.95,
              color='k',edgecolor='none')
for d in deletions_middle:
  plt.axvspan(d[0],d[1],ymin=0.925,ymax=0.95,
              color='white',edgecolor='none')
for m in matches_middle:
  plt.axvspan(m[0],m[1],ymin=0.925,ymax=0.95,
              color=blue,alpha=0.4,edgecolor='none')

#do spanners cter shift
for ins in insertions_cter:
  plt.axvspan(ins[0],ins[1],ymin=0.9,ymax=0.925,
              color='k',edgecolor='none')
for d in deletions_cter:
  plt.axvspan(d[0],d[1],ymin=0.9,ymax=0.925,
              color='white',edgecolor='none')
for m in matches_cter:
  plt.axvspan(m[0],m[1],ymin=0.9,ymax=0.925,
              color=orange,alpha=0.4,edgecolor='none')

#plot the transmembrane helix representation
for sp in membrane_spanners:
  plt.axvspan(sp[0],sp[1],ymin=0.0,ymax=0.9,alpha=0.1,color='k')

#zero line
dots_x = np.linspace(rnums[0],rnums[-1],200)
dots_y = [0] * 200

plt.plot(dots_x,dots_y,'o',color='k',markersize=0.5)

plt.ylim((-0.3,0.3))
plt.xlim((rnums[0],rnums[-1]))
plt.xlabel('Residue Number',fontsize='large')
plt.ylabel('\u0394QMEANBrane Score', fontsize='large')

plt.savefig("alignment_comparison.png",dpi=300)
