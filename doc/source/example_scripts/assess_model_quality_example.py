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

# we now have everything we need and obtain a QMEANScorer object
# to estimate the classic QMEAN4 score
# DISCLAIMER: to use the full power of QMEAN, you should 
# additionally add distance constraints!
scorer = QMEANScorer(crambin, psipred = psipred_handler, 
                     accpro= accpro_handler)

# we get the qmean4 score out and dump the according reference 
# set plot
print("QMEAN4:", scorer.qmean4_score)
scorer.QMEAN4ReferencePlot("qmean_4_ref_plot.png")

# local scores are available as a dictionary
print("Local Scores:", scorer.local_scores)
