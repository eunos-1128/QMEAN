import os
from ost.io import LoadPDB
from qmean import QMEANScorer
from qmean import PSIPREDHandler
from qmean import ACCPROHandler

# we use our favourite example, crambin to demonstrate the usage 
# of the AssessModelQuality function
crambin_path = os.path.join("example_data","1CRN.pdb")
crambin = LoadPDB(crambin_path)

# We define the psipred and accpro prediction in a data dict,
# You can get them for example from:
# PSIPRED: http://bioinf.cs.ucl.ac.uk/psipred/
# ACCPRO: http://scratch.proteomics.ics.uci.edu/

data = dict()

# the sequence
data["seq"] = "TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN"

# the stuff from PSIPRED
data["ss"] = "CCCCCCHHHHHHHCCCCCCCCCHHHHHHCCCCEECCCCCCCCCCCC"
data["conf"] = "9888980234220123269989453431158088299999999989"

# the stuff from ACCPRO
data["acc"] = "ebbbeeeebbeebeebeeeeeeeeebbeebebeebeeeeeeeeeee"

# create the handler
psipred_handler = PSIPREDHandler(data) 
accpro_handler = ACCPROHandler(data)

# we now have everything we need and obtain a QMEANScorer 
# object to estimate the classic QMEAN4 score
# DISCLAIMER: to use the full power of QMEAN, you should 
# additionally add distance constraints!
scorer = QMEANScorer(crambin, psipred = psipred_handler, 
                     accpro= accpro_handler)

# we get the qmean4 and qmean6 score, as well as their 
# Z-score counterparts 
print("QMEAN4:", scorer.qmean4_score)
print("QMEAN6:", scorer.qmean6_score)
print("QMEAN4 Z-score:", scorer.qmean4_z_score)
print("QMEAN6 Z-score:", scorer.qmean6_z_score)

# local scores are available as a dictionary
print("Local Scores:", scorer.local_scores)

# with the rise of QMEANDisCo, average local score has been 
# introduced as a global score predictor. This also comes with 
# an error. However, those errors make only sense when DisCo 
# is provided. We therefore get NaN here...
print("Avg. Local Score:", scorer.avg_local_score)
print("Avg. Local Score Error:", scorer.avg_local_score_error)

# the assessed model is a view of our input structure that we 
# can access as well. If we want to read out local scores from 
# there, we first need to map them.
scorer.AssignModelBFactors()
model = scorer.model
for r in model.residues:
  print(r, r.atoms[0].b_factor)

# the results can also be visualized with various plots
scorer.LocalProfilePlot("local_profile.png")
scorer.QMEAN4ReferencePlot("qmean4_ref_plot.png")
scorer.QMEAN6ReferencePlot("qmean6_ref_plot.png")
scorer.QMEAN4SliderPlot("qmean4_sliders")
scorer.QMEAN6SliderPlot("qmean6_sliders")
