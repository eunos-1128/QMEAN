# Import the required modules
from ost.io import LoadPDB
from ost.seq.alg import AlignToSEQRES
from qmean import AssessModelQuality
from qmean import DCScore
from qmean import LoadDCData

# Let's compute the scores for some models of an example protein (4J32 chain A)
model01 = LoadPDB("example_data/4j32_A_model01.pdb")
model05 = LoadPDB("example_data/4j32_A_model05.pdb")

# Select only the peptide
model01_sele = model01.Select('peptide = true and ligand = false')
model05_sele = model05.Select('peptide = true and ligand = false')

# First the DCData object is loaded
dc_data_4j32_A = LoadDCData("example_data/4j32_A.dat")

# For computation of DisCo scores the residue sequence of the model needs
# to be aligned with the SEQRES and the model has to be attachted to the 
# alignment
model01_aln = AlignToSEQRES(model01_sele, dc_data_4j32_A.seq)
model01_aln.AttachView(1,model01_sele) 

model05_aln = seq.alg.AlignToSEQRES(model05_sele, dc_data_4j32_A.seq)
model05_aln.AttachView(1,model05_sele) 

# Now we can compute the DisCo scores for our two models ...
model01_score = DCScore(model01_aln, dc_data_4j32_A)
model05_score = DCScore(model05_aln, dc_data_4j32_A)

print "DisCo score of model 01:"
print model01_score
print "\nDisCo score of model 05:"
print model05_score

# ... or we can directly compute QMEANDisCo scores
out_path = os.path.join("example_out","qmeandisco_output","model01")
AssessModelQuality(model01_sele, output_dir = out_path, plots = False,
                   global_scores = False, dc = dc_data_4j32_A)

out_path = os.path.join("example_out","qmeandisco_output","model05")
AssessModelQuality(model05_sele, output_dir = out_path, plots = False,
                   global_scores = False, dc = dc_data_4j32_A)
