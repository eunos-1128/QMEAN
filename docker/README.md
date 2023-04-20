A Container for QMEAN Protein Structure Model Evaluation
========================================================

This Dockerfile describes an image of the QMEAN software package ([service](
https://swissmodel.expasy.org/qmean/),
[Git](https://git.scicore.unibas.ch/schwede/QMEAN),
[Docker](https://git.scicore.unibas.ch/schwede/QMEAN/container_registry)).

Table Of Contents
-----------------

* [Available Scoring Functions](#scoringfunctions)

* [Input Requirements](#inputrequirements)

* [Structure Processing](#structureprocessing)

* [Obtain the image (Docker `pull`)](#qmeanpull)

* [Obtain an optimised image (Docker `build`)](#qmeanbuild)

* [Additional requirements](#additionalrequirements)

* [Score](#score)

* [Singularity](#singularity)

* [Results](#results)

* [Examples](#examples)

<a name="scoringfunctions"></a>Available Scoring Functions
----------------------------------------------------------

The following scoring functions are implemented:

[QMEANDisCo](https://doi.org/10.1093/bioinformatics/btz828): 
> Studer, G., Rempfer, C., Waterhouse, A.M., Gumienny, R., Haas, J., Schwede, T. QMEANDisCo-distance constraints applied on model quality estimation, Bioinformatics 36, 1765–1771 (2020).
> 
[QMEAN](https://doi.org/10.1093/bioinformatics/btq662):
> Benkert, P., Biasini, M., Schwede, T. Toward the estimation of the absolute quality of individual protein structure models. Bioinformatics 27, 343-350 (2011).
> 
[QMEANBrane](https://doi.org/10.1093/bioinformatics/btu457):
> Studer, G., Biasini, M., Schwede, T. Assessing the local structural quality of transmembrane protein models using statistical potentials (QMEANBrane), Bioinformatics 30, i505–i511 (2014).
> 

For a short description of the different scoring functions refer to the QMEAN 
server [help page](https://swissmodel.expasy.org/qmean/help).

<a name="inputrequirements"></a>Input Requirements
--------------------------------------------------

The container can read protein structures in PDB (.pdb) or mmCIF (.cif/.mmcif) 
format. Compressed files are accepted, too (.gz suffix, e.g. 1crn.pdb.gz).

Providing SEQRES sequences (via the option `--seqres`) is recommended. SEQRES
describes the actual protein sequence of consecutively connected amino acids.
It is not necessarily completely covered by the structural data, think of
missing terminii in a model or gaps. If not provided, sequences are extracted
from the structural data. SEQRES must be provided as one or multiple sequences
in a single FASTA file. 

* One sequence: All chains of the input structure(s) must align to it, 
                the structure(s) are either monomers or homo-oligmers.  
* Multiple sequences: Required for hetero-oligomers, uses name based 
                     chain/SEQRES matching

The container calculates sequence profiles using 
HHblits ([DOI](https://doi.org/10.1186/s12859-019-3019-7),
[Git](https://github.com/soedinglab/hh-suite)) for each unique SEQRES sequence.
If you already have the respective profiles available in A3M format, you can
speed things up (via the option `--profiles`). This only works if you also
provide SEQRES as an input and the master sequence for each profile must match
one of the SEQRES sequences.

<a name="structureprocessing"></a>Structure Processing
------------------------------------------------------

Structures are processed with the **MOL**ecule **C**hec**K**er 
([molck](https://openstructure.org/docs/dev/mol/alg/molalg/?highlight=molck#ost.mol.alg.Molck))
before scoring. In detail:

* Non-standard residues are mapped to their standard counterparts if possible
  (e.g. Phospho-Tyrosine to Tyrosine, Seleno-Methionine to Methionine, etc.).
  Mapping information is derived from the 
  [component dictionary](http://www.wwpdb.org/data/ccd)
* Hydrogen atoms are stripped
* Atoms with zero occupancy are stripped
* Everything except the 20 standard proteinogenic amino acids is stripped
  (after potential mapping of non-standard residues above)
* Unknown atoms are stripped, i.e. atoms that are not expected based on the
  [component dictionary](http://www.wwpdb.org/data/ccd), as found in some force fields.
* Chains are potentially renamed, i.e. if a chain name is ' ', it gets
  assigned a valid chain name.


<a name="qmeanpull"></a>Obtain the image (Docker `pull`)
--------------------------------------------------------

An already built copy of the image for the current
[Dockerfile](docker/Dockerfile) is available in our [GitLab registry](
https://git.scicore.unibas.ch/schwede/QMEAN/container_registry). It comes
without CPU optimisations, so it should run on most CPUs. The image can be
downloaded by either pulling it yourself or let Docker pull it first time you
run it. To actively pull, use the following command:

```terminal
$ docker pull registry.scicore.unibas.ch/schwede/qmean:latest
latest: Pulling from schwede/qmean
...
Digest: sha256:db53a753d46b2525051478f9fa273df2b47a69100663eb70d538b461d52743d5
Status: Downloaded newer image for registry.scicore.unibas.ch/schwede/qmean:latest
registry.scicore.unibas.ch/schwede/qmean:latest
$
```

<a name="qmeanbuild"></a>Obtain an optimised image (Docker `build`)
--------------------------------------------------------

The available [Docker image](
https://git.scicore.unibas.ch/schwede/QMEAN/container_registry) only uses the
most common CPU extensions to run HHblits ([DOI](
https://doi.org/10.1186/s12859-019-3019-7), [Git](
https://github.com/soedinglab/hh-suite)). Modern CPUs can do better (run
HHblits faster), providing AVX2 extensions. To make use of them, the QMEAN
image needs to be rebuild with the Docker build argument `OPT_FOR_CPU`
activated, on the target host. Execute the following command inside the
directory storing the Dockerfile
for QMEAN:

```terminal
$ docker build --build-arg OPT_FOR_CPU=1 -t qmean:avx2 .
...
=> => writing image sha256:93da5aa6613847fb7f8b71ed635fe99f72e7dd58f32b40e1e73e3429547c9d25            0.0s
 => => naming to qmean:avx2
$
```

<a name="additionalrequirements"></a>Additional requirements
------------------------------------------------------------

We need the non-redundant 
[UniClust30 sequence database](https://uniclust.mmseqs.com/) to build sequence
profiles with HHblits. The following files are required:

* `X_a3m.ffdata`
* `X_a3m.ffindex`
* `X_hhm.ffdata`
* `X_hhm.ffindex`
* `X_cs219.ffdata`
* `X_cs219.ffindex`

with `X` being your UniClust30 version of choice. The productive QMEAN server uses
UniClust30 [August 2018](http://wwwuser.gwdg.de/~compbiol/uniclust/2018_08/).
The directory with the files must be mounted to the container:

```terminal
-v <PATH_TO_LOCAL_UNICLUST>:/uniclust30
```

The QMEANDisCo scoring function additionally requires a template library
(QMEAN Template Library aka QMTL) available [here](
https://swissmodel.expasy.org/repository/download/qmtl/qmtl.tar.bz2). The QMTL
is updated weekly following the PDB update cycle. Checking for new versions on
a Thursday is a good idea. **An up to date QMTL is crucial for optimal QMEANDisCo performance.** 
The following files are required:

* `smtl_uniq_cs219.ffdata`
* `smtl_uniq_cs219.ffindex`
* `smtl_uniq_hhm.ffdata`
* `smtl_uniq_hhm.ffindex`
* `CHAINCLUSTERINDEX`
* `indexer.dat`
* `seqres_data.dat`
* `atomseq_data.dat`
* `ca_pos_data.dat`
* `VERSION`

Again, the corresponding directory must be mounted:

```terminal
-v <PATH_TO_LOCAL_QMTL>:/qmtl
```

<a name="score"></a>Score
-------------------------

Having everything setup, you can score `model.pdb` with SEQRES data stored in
`seqres.fasta` using QMEANDisCo:

```terminal
docker run --workdir $(pwd) -v $(pwd):$(pwd) -v <PATH_TO_LOCAL_UNICLUST>:/uniclust30 -v <PATH_TO_LOCAL_QMTL>:/qmtl registry.scicore.unibas.ch/schwede/qmean:latest run_qmean.py model.pdb --seqres seqres.fasta
```

Additionally to the mounts specified above, the current working directory 
containing the input files is mounted into the container and specified 
as workdir.

The following gives more details on additional command line arguments:

```terminal
docker run registry.scicore.unibas.ch/schwede/qmean:latest run_qmean.py --help
```

<a name="singularity"></a>Singularity
-------------------------------------

A Singularity Image can directly be pulled & build from our registry:

```terminal
singularity build qmean_container.sif docker://registry.scicore.unibas.ch/schwede/qmean:latest
```

Singularity allows to directly access the current working directory from within the container,
so scoring simplifies to:

```terminal
singularity run -B <PATH_TO_LOCAL_UNICLUST>:/uniclust30 -B <PATH_TO_LOCAL_QMTL>:/qmtl qmean_container.sif run_qmean.py model.pdb  --seqres seqres.fasta
```

<a name="results"></a>Results
-----------------------------

Results are JSON formatted. For each model there is an entry with following 
keys:

* `chains`: Contains the ATOMSEQ/SEQRES mapping for each chain. The ATOMSEQ is 
  extracted from the structure, SEQRES is provided by the user. If SEQRES is not 
  provided, SEQRES equals the ATOMSEQ
* `original_name`: Filename of the input model
* `preprocessing`: Summarizes the outcome of input structure processing described 
  above
* `scores`: Model specific scores. No matter what scoring method you're running 
  (QMEAN [1] /QMEANDisCo [2] /QMEANBrane [3]), you have two keys: 
  `local_scores` and `global_scores`. While the data in `global_scores` 
  calculates values according [1], the `local_scores` are method dependent. 
  For `global_scores` you get another dictionary with keys:

  * `acc_agreement_norm_score`
  * `acc_agreement_z_score`
  * `avg_local_score` (mode dependent)
  * `avg_local_score_error` (only set for [2])
  * `cbeta_norm_score`
  * `cbeta_z_score`
  * `nteraction_norm_score`
  * `interaction_z_score`
  * `packing_norm_score`
  * `packing_z_score`
  * `qmean4_norm_score`
  * `qmean4_z_score`
  * `qmean6_norm_score`
  * `qmean6_z_score`
  * `ss_agreement_norm_score`
  * `ss_agreement_z_score`
  * `torsion_norm_score`
  * `torsion_z_score`

  For `local_scores` you get another dictionary with chain names as keys and 
  lists with local scores as value. The local score lists have the same length 
  as the corresponding SEQRES (ATOMSEQ if SEQRES is not given as an input). 
  The location of a residue in SEQRES determines the location of its local score 
  in the score list.
* `qmeanbrane_membrane`: Only available in case of QMEANBrane. Defines
  transmembrane residues with the same SEQRES mapping as the local scores.


<a name="examples"></a>Examples
-------------------------------

Example data to run all available scoring functions are available in the
[qmean_qmeandisco_example](docker/qmean_qmeandisco_example) and
[qmeanbrane_example](docker/qmeanbrane_example) directory.


[comment]: <> ( LocalWords:  QMEANDisCo mmCIF JSON GitLab DBeacons cd OST )
[comment]: <> ( LocalWords:  schwede qmean sha )
