A QMEAN DisCo Container for the 3D Beacon's Project
===================================================

This Dockerfile describes an image of QMEAN ([service](
https://swissmodel.expasy.org/qmean/),
[Git](https://git.scicore.unibas.ch/schwede/QMEAN),
[Docker](https://git.scicore.unibas.ch/schwede/QMEAN/container_registry))
dedicated to the 3D Beacons project ([Git](https://github.com/3D-Beacons)).
More specifically its a containerised QMEAN DisCo instance, helping the
3D Beacons back-end to get common scores for structural models.

Per se, the image and QMEAN script can be used as foundation for different QMEAN
applications. But since a general QMEAN app is out of scope here, the
functionality as-is is limited in terms of input, output and QMEAN scores
calculated.

What this container will do, is simply take a [mmCIF](http://mmcif.rcsb.org)
file, calculate the QMEANDisCo score of the peptides found in the input and
return the local & global scores in JSON format as output.


Obtain the image (Docker `pull`)
------------------------------

An already built copy of the image for the current Dockerfile is available in
our [GitLab registry](
https://git.scicore.unibas.ch/schwede/QMEAN/container_registry). It can be
downloaded by either pulling it yourself or let Docker pull it first time you
run it. To actively pull use the following command:

```terminal
$ docker pull registry.scicore.unibas.ch/schwede/qmean:4.0.0-3d-Beacons
4.0.0-3d-Beacons: Pulling from schwede/qmean
...
Digest: sha256:bedafe6c91a2b51197626640e2db3679b604cd794f2b3ba41bc36bead8727c7d
Status: Downloaded newer image for registry.scicore.unibas.ch/schwede/qmean:3DBeacons
registry.scicore.unibas.ch/schwede/qmean:4.0.0-3d-Beacons
$
```



[comment]: <> ( LocalWords:  QMEAN DisCo mmCIF JSON GitLab DBeacons cd OST )
[comment]: <> ( LocalWords:  schwede qmean sha )
