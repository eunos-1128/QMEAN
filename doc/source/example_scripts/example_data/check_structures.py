aln_files = ["1pvv.1.A.fasta",
             "1v1v.1.A.fasta",
             "2w37.1.A.fasta",
             "3grf.1.A.fasta",
             "4nf2.1.A.fasta"]

tpl_files = ["1pvv.1.A.pdb",
             "1v1v.1.A.pdb",
             "2w37.1.A.pdb",
             "3grf.1.A.pdb",
             "4nf2.1.A.pdb"]


for af, tf in zip(aln_files, tpl_files):
  print "processing",af,tf

  aln = io.LoadAlignment(af)
  tpl = io.LoadPDB(tf)

  tpl_seq = aln.GetSequence(1)

  if len(tpl_seq.GetGaplessString()) != len(tpl.residues):
    raise RuntimeError("AAAAAAAAAAAAAAAA")

  for a,b in zip(tpl_seq.GetGaplessString(), tpl.residues):
    if a != b.one_letter_code:
      raise RuntimeError("AAAAAAAAAA")

