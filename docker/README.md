A QMEANDisCo Container for the 3D Beacon's Project
===================================================

This Dockerfile describes an image of QMEAN ([service](
https://swissmodel.expasy.org/qmean/),
[Git](https://git.scicore.unibas.ch/schwede/QMEAN),
[Docker](https://git.scicore.unibas.ch/schwede/QMEAN/container_registry))
dedicated to the 3D Beacons project ([Git](https://github.com/3D-Beacons)).
More specifically its a containerised [QMEANDisCo](
https://doi.org/10.1093/bioinformatics/btz828) instance, helping the 3D Beacons
back-end to get common scores for structural models.

Per se, the image and QMEAN script can be used as foundation for different QMEAN
applications. But since a general QMEAN app is out of scope here, the
functionality as-is is limited in terms of input, output and QMEAN scores
calculated.

What this container will do, is simply take a [mmCIF](http://mmcif.rcsb.org)
file, calculate the [QMEANDisCo](
https://doi.org/10.1093/bioinformatics/btz828) score of the peptides found in the
input and return the local & global scores in JSON format as output.


<a name="qmeanpull"></a>Obtain the image (Docker `pull`)
------------------------------

An already built copy of the image for the current Dockerfile is available in
our [GitLab registry](
https://git.scicore.unibas.ch/schwede/QMEAN/container_registry). It can be
downloaded by either pulling it yourself or let Docker pull it first time you
run it. To actively pull use the following command:

```terminal
$ docker pull registry.scicore.unibas.ch/schwede/qmean:4.2.0-3d-Beacons
4.2.0-3d-Beacons: Pulling from schwede/qmean
...
Digest: sha256:bedafe6c91a2b51197626640e2db3679b604cd794f2b3ba41bc36bead8727c7d
Status: Downloaded newer image for registry.scicore.unibas.ch/schwede/qmean:3DBeacons
registry.scicore.unibas.ch/schwede/qmean:4.2.0-3d-Beacons
$
```


Obtain the image (Docker `build`)
------------------------------

In case you have to build your own container, e.g. to adjust the user ID running
QMEAN inside the container to match file ownership, all you need is the
[Dockerfile](docker/Dockerfile) from the QMEAN Git repository.

If you are in the same directory as the Dockerfile, you can build the container
like this:

```terminal
$ docker build -t registry.scicore.unibas.ch/schwede/qmean:4.2.0-3d-Beacons .
...
$
```

Switching the user ID/ group ID of the user inside the container, use build
arguments:

```terminal
$ docker build --build-arg QMEAN_USER_ID=$(id -u) --build-arg QMEAN_GROUP_ID=$(id -g) -t registry.scicore.unibas.ch/schwede/qmean:4.2.0-3d-Beacons .
...
$
```

That command will use the IDs of the current user who builds the container.
Since the user inside the container is created at the end of the Dockerfile,
[pulling](#qmeanpull) the container before the build from the QMEAN GitLab
registry will make the build almost instantly.

Other build arguments available are:

- `CPUS_FOR_MAKE`, processes to be started by make, parameter to `-j`

- `OST_VERSION` version (Git tag) of OpenStructure ([website](
  https://openstructure.org), [Git](
  https://git.scicore.unibas.ch/schwede/openstructure)) to be used in the build

- `QMEAN_VERSION` version (Git tag) of QMEAN to be build


[comment]: <> ( LocalWords:  QMEANDisCo mmCIF JSON GitLab DBeacons cd OST )
[comment]: <> ( LocalWords:  schwede qmean sha )
