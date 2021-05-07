Run QMEAN/QMEANDisCo
====================

The simplest you can do is simply run QMEAN on the provided model:

```terminal
sudo docker run --workdir $(pwd) -v $(pwd):$(pwd) -v <PATH_TO_LOCAL_UNICLUST>:/uniclust30 registry.scicore.unibas.ch/schwede/qmean:4.2.0 run_qmean.py model.pdb --method QMEAN
```

However, if possible, you should add the SEQRES sequence which is available in 
targets.fasta:

```terminal
sudo docker run --workdir $(pwd) -v $(pwd):$(pwd) -v <PATH_TO_LOCAL_UNICLUST>:/uniclust30 registry.scicore.unibas.ch/schwede/qmean:4.2.0 run_qmean.py model.pdb --method QMEAN --seqres targets.fasta
```
The directory contains pre-computed profiles that match the 
sequences in targets.fasta. To speed things up, you can provide them as argument:

```terminal
sudo docker run --workdir $(pwd) -v $(pwd):$(pwd) -v <PATH_TO_LOCAL_UNICLUST>:/uniclust30 registry.scicore.unibas.ch/schwede/qmean:4.2.0 run_qmean.py model.pdb --method QMEAN --seqres targets.fasta --profiles query_hhblits_one.a3m query_hhblits_two.a3m
```

Normally, the container creates a temporary directory to store intermediate 
results. If you want to investigate them, or access the sequence profiles,
you can provide a working dir as argument (must be absolute path). 
A directory gets created in the working dir for each unique SEQRES named after 
the md5 hash of the SEQRES:

```terminal
sudo docker run --workdir $(pwd) -v $(pwd):$(pwd) -v <PATH_TO_LOCAL_UNICLUST>:/uniclust30 registry.scicore.unibas.ch/schwede/qmean:4.2.0 run_qmean.py model.pdb --method QMEAN --seqres targets.fasta --profiles query_hhblits_one.a3m query_hhblits_two.a3m --workdir $(pwd)/my_workdir
```

So far we only computed the QMEAN scoring function. QMEANDisCo includes 
structural information and significantly improves prediction performance 
of per-residue scores. QMEANDisCo is the default value of the --method
argument. We nevertheless explicitely specify it here. QMEANDisCo adds
the additonal requirement of the QMTL (QMEAN Template Library) which needs
to be mounted.

```terminal
sudo docker run --workdir $(pwd) -v $(pwd):$(pwd) -v <PATH_TO_LOCAL_UNICLUST>:/uniclust30 -v <PATH_TO_LOCAL_QMTL>:/qmtl registry.scicore.unibas.ch/schwede/qmean:4.2.0 run_qmean.py model.pdb --method QMEANDisCo --seqres targets.fasta --profiles query_hhblits_one.a3m query_hhblits_two.a3m --workdir my_workdir
```

