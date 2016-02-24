# Import the required modules
from qmean.dc_data import TemplateInformation
from qmean import FillDCData
from ost.io import LoadAlignment
from ost.io import LoadSequence
from ost.seq.alg import BLOSUM62


# In this example we compute the DCData object for the SEQRES of 4J32 (chain A)
# a protein from Salmonella Typhimurium

# Load the SEQRES and convert it to a string
seqres = LoadSequence("example_data/4j32_A_seqres.fasta")
seqres_as_string = seqres.GetGaplessString()

# Load a multiple sequence alignment of the SEQURES with homologous
# template sequences
ms_aln = LoadAlignment("example_data/4j32_A_msaln.fasta")

print ms_aln

# Provide a list of paths to template structures. Paths need to be
# in the same order as the templates are in the multiple sequence
# alignment
templates = ["example_data/4j32_A_tpls/2k8f_B.pdb.gz",
             "example_data/4j32_A_tpls/4bi3_B.pdb.gz",
             "example_data/4j32_A_tpls/4bi4_A_2.pdb.gz",
             "example_data/4j32_A_tpls/4hfk_A.pdb.gz",
             "example_data/4j32_A_tpls/4bi3_A.pdb.gz",
             "example_data/4j32_A_tpls/4bi4_A_1.pdb.gz",
             "example_data/4j32_A_tpls/4bi8_B.pdb.gz",
             "example_data/4j32_A_tpls/4hfl_A.pdb.gz"]

# Cluster templates by sequence similarity
tpl_inf = TemplateInformation(seqres_as_string)
cluster_list = tpl_inf.ClusterTemplates(ms_aln)

# Extract the information from homologous templates and generate
# the DCData object
disco_data = FillDCData(ms_aln, cluster_list, templates, BLOSUM62)

# Let's look at some of the information  
print "SEQRES of 4J32 chain A:" 
print disco_data.seq
print "\nAverage sequence similarities to SEQRES of templates in same clusters"
print disco_data.cluster_seq_sim
print "\nAverage sequence identities to SEQRES of templates in same clusters"
print disco_data.cluster_seq_ids







