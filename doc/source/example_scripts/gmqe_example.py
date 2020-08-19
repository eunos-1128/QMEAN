from ost import io
from ost import seq
from ost import mol
from qmean import PSIPREDHandler
from qmean import DisCoContainer
from qmean import GMQE

# Let's model crambin! 
trg_sequence = seq.CreateSequence('crambin', 'TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN')
trg_profile = io.LoadSequenceProfile('example_data/1crn.1.A.hhm')
dc = DisCoContainer.Load('example_data/1crn_dc.dat')

# Setup PSIPREDHandler
data = dict()
data["seq"] = str(trg_sequence)
data["ss"] = 'CCCCCCHHHHHHHCCCCCCCCCHHHHHHCCCCEECCCCCCCCCCCC'
data["conf"] = '9888980234220123269989453431158088299999999989'
psipred_handler = PSIPREDHandler(data) 

# GMQE is calculated by using the structure with PDB ID 3szs as template
seqres_aln = io.LoadAlignment('example_data/1crn_3szs_aln.fasta')
tpl = io.LoadPDB('example_data/3szs.2.A.pdb')
tpl_profile = io.LoadSequenceProfile('example_data/3szs.2.A.hhm')

# Assign secondary structure, the ss_agreement score is invalid otherwise
mol.alg.AssignSecStruct(tpl)

# The template has no gaps, the atomseq alignment therefore matches the
# seqres alignment. Be aware that the AttachView function doesn't check
# whether sequence and attached view match. 
aln = seqres_aln
aln.AttachView(1, tpl.CreateFullView())

# Generate target specific GMQE object and estimate GMQE
# For optimal performance you need to provide a context profile database
# (see documentation)
scorer = GMQE(trg_sequence, psipred_handler, disco = dc, profile = trg_profile)
gmqe_result = scorer.PredictGMQE(aln, seqres_aln, tpl_profile = tpl_profile)
print('Expected model quality when using 3szs as template: ', gmqe_result[0])
print('Scores: ', gmqe_result[1])

