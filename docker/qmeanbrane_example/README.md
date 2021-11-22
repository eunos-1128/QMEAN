Run QMEANBrane
==============

This example runs QMEANBrane on the models of the example application in the 
[QMEANBrane manuscript](https://doi.org/10.1093/bioinformatics/btu457).

We don't use QMEANDisCo, so the QMTL is not required. The most basic docker 
command that generates output in out.json would be:

```terminal
sudo docker run --workdir $(pwd) -v $(pwd):$(pwd) -v <PATH_TO_LOCAL_UNICLUST>:/uniclust30 registry.scicore.unibas.ch/schwede/qmean:4.2.0 run_qmean.py original_hhblits_alignment.pdb shift_in_front_helix_four.pdb shift_into_middle.pdb shift_towards_cter.pdb --method QMEANBrane
```

However, if possible, you should add the SEQRES sequence which is available in 
targets.fasta:

```terminal
sudo docker run --workdir $(pwd) -v $(pwd):$(pwd) -v <PATH_TO_LOCAL_UNICLUST>:/uniclust30 registry.scicore.unibas.ch/schwede/qmean:4.2.0 run_qmean.py original_hhblits_alignment.pdb shift_in_front_helix_four.pdb shift_into_middle.pdb shift_towards_cter.pdb --method QMEANBrane --seqres targets.fasta
```
The directory contains an already pre-computed profile which matches the 
sequence in targets.fasta. To speed things up, you can provide it as argument:

```terminal
sudo docker run --workdir $(pwd) -v $(pwd):$(pwd) -v <PATH_TO_LOCAL_UNICLUST>:/uniclust30 registry.scicore.unibas.ch/schwede/qmean:4.2.0 run_qmean.py original_hhblits_alignment.pdb shift_in_front_helix_four.pdb shift_into_middle.pdb shift_towards_cter.pdb --method QMEANBrane  --seqres targets.fasta --profiles query_hhblits.a3m
```

Normally, the container creates a temporary directory to store intermediate 
results. If you want to investigate them, or access the sequence profiles,
you can provide a working dir as argument (must be absolute path). 
A directory gets created in the working dir for each unique SEQRES named 
after the md5 hash of the SEQRES:

```terminal
sudo docker run --workdir $(pwd) -v $(pwd):$(pwd) -v <PATH_TO_LOCAL_UNICLUST>:/uniclust30 registry.scicore.unibas.ch/schwede/qmean:4.2.0 run_qmean.py original_hhblits_alignment.pdb shift_in_front_helix_four.pdb shift_into_middle.pdb shift_towards_cter.pdb --method QMEANBrane  --seqres targets.fasta --profiles query_hhblits.a3m --workdir $(pwd)/my_workdir
```

