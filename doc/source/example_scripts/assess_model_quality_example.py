#import the required modules
from qmean import AssessModelQuality
from qmean import PSIPREDHandler
from qmean import ACCPROHandler
from ost.io import LoadPDB
import os

#we use our favourite example, crambin to demonstrate the usage of
#the AssessModelQuality function
crambin_path = os.path.join("example_data","1CRN.pdb")
crambin = LoadPDB(crambin_path)

#We define the psipred and accpro prediction in a data dict,
#You can get them for example from:
#PSIPRED: http://bioinf.cs.ucl.ac.uk/psipred/
#ACCPRO: http://scratch.proteomics.ics.uci.edu/

data = dict()

#the sequence
data["seq"] = "TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN"

#the stuff from PSIPRED
data["ss"] = "CCCCCCHHHHHHHCCCCCCCCCHHHHHHCCCCEECCCCCCCCCCCC"
data["conf"] = "9888980234220123269989453431158088299999999989"

#the stuff from ACCPRO
data["acc"] = "ebbbeeeebbeebeebeeeeeeeeebbeebebeebeeeeeeeeeee"

#create the handler
psipred_handler = PSIPREDHandler(data) 
accpro_handler = ACCPROHandler(data)

#we now have everything we need and can run the quality assessment
out_path = os.path.join("example_out","crambin_assessment")
AssessModelQuality(crambin, output_dir=out_path,
                   psipred = psipred_handler, accpro= accpro_handler)
