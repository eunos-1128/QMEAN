from qmean import DisCoContainer
import numpy as np
import matplotlib.pyplot as plt
import os

aln_files = ["example_data/1pvv.1.A.fasta", 
             "example_data/1v1v.1.A.fasta", 
             "example_data/2w37.1.A.fasta", 
             "example_data/3grf.1.A.fasta",
             "example_data/4nf2.1.A.fasta"]

tpl_files = ["example_data/1pvv.1.A.pdb", 
             "example_data/1v1v.1.A.pdb",
             "example_data/2w37.1.A.pdb", 
             "example_data/3grf.1.A.pdb",
             "example_data/4nf2.1.A.pdb"]

alns = list()
tpls = list()

for af, tf in zip(aln_files, tpl_files):
  alns.append(io.LoadAlignment(af))
  tpls.append(io.LoadPDB(tf))

# First sequence of all alns is the same: the SEQRES of the model
seqres = io.LoadSequence("example_data/target_seq.fasta")

dc = DisCoContainer(seqres)

for aln, tpl in zip(alns, tpls):

  # attach the template to the alignment
  # there is NO CHECK whether the template sequence and the 
  # structure match, you have to make sure of that by yourself. 
  # The data for this example is preprocessed accordingly.
  aln.AttachView(1, tpl.CreateFullView())

  # in here we store the positions
  positions = geom.Vec3List()

  # in here the mapping to the seqres
  mapping = list()

  current_rnum = 0
  for col in aln:

    if col[0] != '-':
      current_rnum += 1

    # check for aligned col
    if col[0] != '-' and col[1] != '-':
      res = col.GetResidue(1)
      ca = res.FindAtom("CA")
      if ca.IsValid():
        positions.append(ca.GetPos())
        mapping.append(current_rnum)

  dc.AddData(aln, positions, mapping)

# after calling that function, the dc object is ready to be used
dc.CalculateConstraints()

# Let's first score a model. The residue numbers in the model
# must match the position in the SEQRES (indexing starting at one)
model = io.LoadPDB("example_data/model_01.pdb")
scores = dc.GetScores(model.CreateFullView())

# the length of the returned list corresponds to the number of
# residues in the model...
print "n residues: ", len(model.residues)
print "n scores:", len(scores)

# we assign the DisCo scores to the bfactors and dump the model
for res, score in zip(model.residues, scores):
  print res, "DisCo score:", score
  for a in res.atoms:
    a.b_factor = score

io.SavePDB(model, "example_out/assigned_disco.pdb")

# we just plot one of the constraints:
rnum_i = 58
rnum_j = 261
c = dc.GetConstraint(rnum_i, rnum_j)
x = np.linspace(0.0, len(c) * dc.GetBinSize(), len(c))
plt.plot(x,c)
plt.xlabel("dist 58-261")
plt.ylabel("score")
plt.savefig(os.path.join("example_out", "distance_constraint.png"))

