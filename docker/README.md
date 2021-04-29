A Container for QMEAN Protein Structure Model Evaluation
========================================================

This Dockerfile describes an image of the QMEAN software package ([service](
https://swissmodel.expasy.org/qmean/),
[Git](https://git.scicore.unibas.ch/schwede/QMEAN),
[Docker](https://git.scicore.unibas.ch/schwede/QMEAN/container_registry))

The following scoring functions are implemented:

* [QMEANDisCo](https://doi.org/10.1093/bioinformatics/btz828)
* [QMEAN](https://doi.org/10.1093/bioinformatics/btq662)
* [QMEANBrane](https://doi.org/10.1093/bioinformatics/btu457)

For a short description of the different scoring functions we refer to the QMEAN 
server [help page](https://swissmodel.expasy.org/qmean/help).

Input Requirements
------------------

The container can read protein structures in PDB (.pdb) or mmCIF (.cif/.mmcif) 
format. Compressed files are accepted too (.gz suffix, i.e. 1crn.pdb.gz).

Providing SEQRES sequences is recommended, this is the actual protein sequence
which is not necessarily fully covered by the structural data. If not provided,
sequences are extracted from the structural data. SEQRES must be provided as one
or several sequences in a single FASTA file. 

* One sequence: All chains of the input structure(s) must align to it, 
                the structures are either monomers or homo-oligmers.  
* Several sequences: Required for hetero-oligomers, uses name based 
                     chain/SEQRES matching

The container calculates sequence profiles using 
[HHblits](https://github.com/soedinglab/hh-suite) for each unique SEQRES 
sequence. If you already have the respective profiles available in a3m format, 
you can speed things up. This only works if you also provide SEQRES as an input
and the master sequence for each profile must match one of the SEQRES sequences.

<a name="qmeanpull"></a>Obtain the image (Docker `pull`)
--------------------------------------------------------

An already built copy of the image for the current Dockerfile is available in
our [GitLab registry](
https://git.scicore.unibas.ch/schwede/QMEAN/container_registry). It can be
downloaded by either pulling it yourself or let Docker pull it first time you
run it. To actively pull use the following command:

```terminal
$ docker pull registry.scicore.unibas.ch/schwede/qmean:4.2.0
4.2.0: Pulling from schwede/qmean
...
Digest: sha256:db53a753d46b2525051478f9fa273df2b47a69100663eb70d538b461d52743d5
Status: Downloaded newer image for registry.scicore.unibas.ch/schwede/qmean:4.2.0
registry.scicore.unibas.ch/schwede/qmean:4.2.0
$
```

Additional requirements 
-----------------------

We need the non-redundant 
[UniClust30 sequence database](https://uniclust.mmseqs.com/) to build sequence
profiles with HHblits. The following_files are required:

* <X>_a3m.ffdata
* <X>_a3m.ffindex
* <X>_hhm.ffdata
* <X>_hhm.ffindex
* <X>_cs219.ffdata
* <X>_cs219.ffindex

with X being your UniClust30 version of choice. The productive QMEAN server uses
UniClust30 [August 2018](http://wwwuser.gwdg.de/~compbiol/uniclust/2018_08/).
The directory with the files must be mounted to the container:

```terminal
-v <PATH_TO_LOCAL_UNICLUST>:/uniclust30
```

The QMEANDisCo scoring function additionally requires a WEEKLY UPDATED template 
library (QMEAN Template Library aka QMTL) available HERE. 
The following files are required:

* smtl_uniq_cs219.ffdata
* smtl_uniq_cs219.ffindex
* smtl_uniq_hhm.ffdata
* smtl_uniq_hhm.ffindex
* CHAINCLUSTERINDEX
* indexer.dat
* seqres_data.dat
* atomseq_data.dat
* ca_pos_data.dat
* VERSION

Again, the according directory must be mounted:

```terminal
-v <PATH_TO_LOCAL_QMTL>:/qmtl
```

Score
-----

Having everything setup, you can score model.pdb with SEQRES data stored in
seqres.fasta using QMEANDisCo:

```terminal
docker run --workdir $(pwd) -v $(pwd):$(pwd) -v <PATH_TO_LOCAL_UNICLUST>:/uniclust30 -v <PATH_TO_LOCAL_QMTL>:/qmtl registry.scicore.unibas.ch/schwede/qmean:4.2.0 qmeandisco model.pdb --seqres seqres.fasta
```

Additionally to the mounts specified above, the current working directory 
containing the input files is mounted into the container and specified 
as workdir.

The following gives more details on additional command line arguments:

```terminal
docker run registry.scicore.unibas.ch/schwede/qmean:4.2.0 qmeandisco --help
```

Singularity
-----------
TODO


[comment]: <> ( LocalWords:  QMEANDisCo mmCIF JSON GitLab DBeacons cd OST )
[comment]: <> ( LocalWords:  schwede qmean sha )
