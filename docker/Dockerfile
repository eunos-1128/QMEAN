FROM centos:7.4.1708

# Arguments, optional
ARG CPUS_FOR_MAKE=1

# Environment variables
ENV CFLAGS="-march=x86-64 -mtune=generic -isystem /usr/local/include/boost/ -isystem /usr/local/openmm/include" \
    HHLIB="/usr/local" \
    LD_LIBRARY_PATH="/usr/local/lib64:/usr/local/lib:$LD_LIBRARY_PATH" \
    PERL5LIB="/usr/local/scripts" \
    PYTHONPATH="/usr/local/lib64/python3.6/site-packages/:${PYTHONPATH}" \
    SCCRE_GTLB="https://git.scicore.unibas.ch/schwede" \
    SRC_FLDR="/tmp" \
    TEMP="/tmp" \
    TEMPDIR="/tmp" \
    TMP="/tmp" \
    TMPDIR="/tmp"

# Define environment variables depending on other variables
ENV CXXFLAGS=${CFLAGS}

# Set up TMP env var. For some reason, it can't be established using ENV in a
# CentOS image.
RUN echo "# For bash users with uid higher" > /etc/profile.d/tmpdir.sh && \
    echo "# than 500 TMPDIR env var points to" >> /etc/profile.d/tmpdir.sh && \
    echo "# /scratch. If $TMPDIR is not" >> /etc/profile.d/tmpdir.sh && \
    echo "# defined." >> /etc/profile.d/tmpdir.sh && \
    echo "# (DISABLED UNTIL WE ACTUALLY " >> /etc/profile.d/tmpdir.sh && \
    echo "# WANT THE SPLIT)" >> /etc/profile.d/tmpdir.sh && \
    echo "#if [ "$UID" -ge 500 ]; then" >> /etc/profile.d/tmpdir.sh && \
    echo "	if [ -z "$TMPDIR" ]; then" >> /etc/profile.d/tmpdir.sh && \
    echo "		export TMPDIR=/tmp" >> /etc/profile.d/tmpdir.sh && \
    echo "	fi" >> /etc/profile.d/tmpdir.sh && \
    echo "	export TMP=${TMPDIR}" >> /etc/profile.d/tmpdir.sh && \
    echo "	export TEMP=${TMPDIR}" >> /etc/profile.d/tmpdir.sh && \
    echo "	export TEMPDIR=${TMPDIR}" >> /etc/profile.d/tmpdir.sh && \
    echo "#fi" >> /etc/profile.d/tmpdir.sh

WORKDIR ${SRC_FLDR}
#
# Compile & install software dependencies
#
# Versions of software used, prefixed with VERSION so $VERSION_<TAB><TAB>
# in the shell produces a list of possible software versions ;)
ENV VERSION_BLAST="2.2.26" \
    VERSION_BOOST="1.68.0" \
    VERSION_CMAKE="3.17.3" \
    VERSION_EIGEN="3.3.1" \
    VERSION_HHBLITS="3.3.0" \
    VERSION_PSIPRED="3.5" \
    VERSION_SSPRO="4.03"
RUN set -eo pipefail; \
    # \
    # Install dependencies for building software \
    # \
    /usr/bin/yum makecache fast; \
    /usr/bin/yum -y install epel-release; \
    /usr/bin/yum update --assumeyes; \
    /usr/bin/yum install --assumeyes bzip2 gcc-c++ compat-libstdc++-33.i686 \
                                     glibc.i686 make openssl-devel pcre-devel \
                                     perl python3 python3-devel zlib-static; \
    /usr/bin/yum --assumeyes clean all; \
    /usr/bin/rm -rf /var/cache/yum; \
    # \
    # Install Python modules \
    # \
    /usr/bin/python3 -m pip install --upgrade --no-cache-dir pip; \
    /usr/bin/python3 -m pip install --no-cache-dir numpy==1.17.4 \
                                                   scipy==1.4.0 \
                                                   matplotlib==3.1.2; \
    # configure matplotlib \
    /usr/bin/mkdir -p $HOME/.config/matplotlib/; \
    echo "backend      : Agg" > $HOME/.config/matplotlib/matplotlibrc; \
    # Import matplot once to generate the fontcache \
    echo "import matplotlib.pyplot" | /usr/bin/python3; \
    # \
    # Install: a recent version of CMake \
    # \
    cd ${SRC_FLDR}; \
    /usr/bin/curl -o cmake.tar.gz -L https://github.com/Kitware/CMake/releases/download/v${VERSION_CMAKE}/cmake-${VERSION_CMAKE}.tar.gz; \
    /usr/bin/mkdir -p cmake/build; \
    /usr/bin/tar --strip-components=1 -C cmake -xf cmake.tar.gz; \
    /usr/bin/rm -f cmake.tar.gz; \
    cd cmake/build; \
    ../bootstrap; \
    /usr/bin/make -j ${CPUS_FOR_MAKE}; \
    /usr/bin/make install; \
    cd ${SRC_FLDR}; \
    /usr/bin/rm -rf cmake; \
    # \
    # Install: Eigen \
    # \
    cd ${SRC_FLDR}; \
    /usr/bin/curl -o eigen.tar.bz2 -L https://gitlab.com/libeigen/eigen/-/archive/${VERSION_EIGEN}/eigen-${VERSION_EIGEN}.tar.bz2; \
    /usr/bin/mkdir eigen/; \
    /usr/bin/tar --strip-components=1 -C eigen/ -xf eigen.tar.bz2; \
    /usr/bin/rm -rf eigen.tar.bz2; \
    /usr/bin/mv eigen/Eigen /usr/local/include; \
    /usr/bin/rm -rf eigen; \
    # \
    # Install: Boost \
    # \
    cd ${SRC_FLDR}; \
    /usr/bin/curl -o boost.tar.bz2 -L https://sourceforge.net/projects/boost/files/boost/${VERSION_BOOST}/boost_${VERSION_BOOST//\./_}.tar.bz2/download; \
    /usr/bin/mkdir boost; \
    /usr/bin/tar --strip-components=1 -C boost -xf boost.tar.bz2; \
    /usr/bin/rm -f boost.tar.bz2; \
    cd boost; \
    # Python 3.6 in CentOS uses pymalloc, hence there is a 'm' at the end of
    # the include dir which Boost 1.68.0 does not consider.
    /usr/bin/sed -i 's/includes ?= \$(prefix)\/include\/python\$(version)/includes ?= \$(prefix)\/include\/python\$(version)m/g' \
                 tools/build/src/tools/python.jam; \
    ./bootstrap.sh --prefix=/usr/local --with-python=/usr/bin/python3; \
    ./b2 --without-mpi; \
    ./b2 install; \
    cd ${SRC_FLDR}; \
    /usr/bin/rm -rf boost; \
    # \
    # Install PsiPred \
    # \
    cd ${SRC_FLDR}; \
    /usr/bin/mkdir psipred; \
    /usr/bin/curl -o psipred.tar.gz -L http://bioinfadmin.cs.ucl.ac.uk/downloads/psipred/old_versions/psipred${VERSION_PSIPRED}.tar.gz; \
    /usr/bin/tar --strip-components=1 -C psipred -xf psipred.tar.gz; \
    /usr/bin/rm -f psipred.tar.gz; \
    cd psipred/src; \
    /usr/bin/make -j ${CPUS_FOR_MAKE}; \
    /usr/bin/cp {chkparse,pfilt,psipass2,psipred,seq2mtx} /usr/local/bin/; \
    cd ..; \
    /usr/bin/cp {runpsipred_single,runpsipred} /usr/local/bin/; \
    /usr/bin/cp -r data /usr/local/; \
    cd ${SRC_FLDR}; \
    /usr/bin/rm -rf psipred; \
    # \
    # Install BLAST 2.2.26 (this is the old version, not BLAST+)\
    # \
    cd ${SRC_FLDR}; \
    /usr/bin/mkdir blast; \
    /usr/bin/curl -o blast.tar.gz -L https://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED/${VERSION_BLAST}/blast-${VERSION_BLAST}-x64-linux.tar.gz; \
    /usr/bin/tar --strip-components=1 -C blast -xf blast.tar.gz; \
    /usr/bin/rm -f blast.tar.gz; \
    /usr/bin/cp blast/bin/* /usr/local/bin/; \
    /usr/bin/mkdir -p /usr/local/data/; \
    /usr/bin/cp blast/data/* /usr/local/data/; \
    /usr/bin/rm -rf blast; \
    # \
    # Install HHblits \
    # \
    cd ${SRC_FLDR}; \
    /usr/bin/mkdir hh-suite; \
    /usr/bin/curl -o hh-suite.tar.gz -L https://github.com/soedinglab/hh-suite/archive/refs/tags/v${VERSION_HHBLITS}.tar.gz; \
    /usr/bin/tar --strip-components=1 -C hh-suite -xf hh-suite.tar.gz; \
    /usr/bin/rm -f hh-suite.tar.gz; \
    cd hh-suite; \
    # Fix paths settings in a Perl module of the HH-Suite
    /usr/bin/sed -i 's/\/cluster\/toolkit\/production\/bioprogs\/psipred\/bin/\/usr\/local\/bin/g' scripts/HHPaths.pm; \
    /usr/bin/sed -i 's/\/cluster\/toolkit\/production\/bioprogs\/psipred\/data/\/usr\/local\/data/g' scripts/HHPaths.pm; \
    /usr/bin/sed -i 's/\/cluster\/toolkit\/production\/bioprogs\/blast\/bin/\/usr\/local\/bin/g' scripts/HHPaths.pm; \
    #/usr/bin/sed -i 's/\/cluster\/databases\/pdb\/all/\/scicore\/data\/databases\/PDB/g' scripts/HHPaths.pm; \
    #/usr/bin/sed -i 's/\/cluster\/databases\/dssp\/bin\/dsspcmbi//g' scripts/HHPaths.pm; \
    /usr/local/bin/cmake . -DCMAKE_INSTALL_PREFIX=/usr/local; \
    /usr/bin/make -j ${CPUS_FOR_MAKE} NO_PNG=1; \
    /usr/bin/make install INSTALL_DIR=/usr/local; \
    cd ${SRC_FLDR}; \
    /usr/bin/rm -rf hh-suite; \
    # \
    # Install SSpro and clean up outdated/ unused files \
    # \
    /usr/bin/mkdir /usr/local/sspro; \
    /usr/bin/curl -o sspro.tar.gz -L http://download.igb.uci.edu/sspro4.tar.gz; \
    /usr/bin/tar --strip-components=1 -C /usr/local/sspro -xf sspro.tar.gz; \
    /usr/bin/rm -f sspro4.tar.gz; \
    cd /usr/local/sspro; \
    /usr/bin/sed -i 's/\/home\/baldig\/jianlin\/sspro4\//\/usr\/local\/sspro\//g' \
                 configure.pl; \
    /usr/bin/perl configure.pl; \
    /usr/bin/rm -rf /usr/local/sspro/data/nr /usr/local/sspro/data/big \
                    /usr/local/sspro/data/pdb_small \
                    /usr/local/sspro/data/pdb_large/old/ \
                    /usr/local/sspro/data/pdb_large/new/ \
                    /usr/local/sspro/blast2.2.8 \
                    /usr/local/sspro/model/sspro \
                    /usr/local/sspro/test \
                    /usr/local/sspro/configure.pl~; \
    /usr/bin/chmod -R og+rx /usr/local/sspro/; \
    # \
    # Cleanup packages only needed to compile/ install things to keep the \
    # image size low. \
    # \
    /usr/bin/yum --assumeyes remove bzip2 gcc gcc-c++ glibc-devel make \
                                    pcre-devel openssl-devel python3-devel \
                                    python3-pip zlib-devel zlib-static; \
    /usr/bin/yum --assumeyes autoremove; \
    /usr/bin/yum --assumeyes clean all; \
    /usr/bin/rm -rf /var/cache/yum; \
    /usr/bin/rm -rf /tmp/*; \
    # removing python3-devel also removes the entire Python 3 installation...
    /usr/bin/yum makecache fast; \
    /usr/bin/yum -y install epel-release; \
    /usr/bin/yum install --assumeyes python3; \
    /usr/bin/yum --assumeyes clean all; \
    /usr/bin/rm -rf /var/cache/yum

#
# Install OpenStructure
#
# Version can be switched via --build-arg OST_VERSION="<TAG|BRANCH>"
ARG OST_VERSION="2.2.0"
ENV VERSION_OPENSTRUCTURE=$OST_VERSION
RUN set -eo pipefail; \
    # \
    # Install dependencies for building software \
    # \
    /usr/bin/yum makecache fast; \
    /usr/bin/yum -y install epel-release; \
    /usr/bin/yum update --assumeyes; \
    /usr/bin/yum install --assumeyes fftw fftw-devel gcc-c++ libpng \
                                          libpng-devel libtiff-devel make \
                                          python3-devel \
                                          sqlite-devel \
    /usr/bin/yum --assumeyes clean all;      \
    /usr/bin/rm -rf /var/cache/yum; \
    cd ${SRC_FLDR}; \
    /usr/bin/mkdir ost; \
    /usr/bin/curl -o ost.tar.gz -L ${SCCRE_GTLB}/openstructure/repository/${OST_VERSION}/archive.tar.gz; \
    /usr/bin/tar --strip-components=1 -C ost -xf ost.tar.gz; \
    /usr/bin/rm -f ost.tar.gz; \
    /usr/bin/mkdir ost/build; \
    cd ost/build; \
    /usr/local/bin/cmake ..  -DCOMPILE_TMTOOLS=OFF \
                             -DENABLE_INFO=OFF \
                             -DENABLE_GFX=OFF \
                             -DENABLE_GUI=OFF \
                             -DUSE_NUMPY=OFF \
                             -DUSE_RPATH=1 \
                             -DOPTIMIZE=1; \
    # 1st run: get everything to build a compounds library
    /usr/bin/make -j ${CPUS_FOR_MAKE} chemdict_tool; \
    # get the compounds dictionary from RCSB PDB
    /usr/bin/curl -o components.cif.gz -L ftp://ftp.wwpdb.org/pub/pdb/data/monomers/components.cif.gz; \
    stage/bin/chemdict_tool create components.cif.gz compounds.chemlib pdb; \
    stage/bin/chemdict_tool update modules/conop/data/charmm.cif \
                                   compounds.chemlib charmm; \
    # 2nd run: build OST with a compounds library
    /usr/local/bin/cmake .. -DCOMPOUND_LIB=compounds.chemlib; \
    /usr/bin/make -j ${CPUS_FOR_MAKE}; \
    /usr/bin/make check; \
    /usr/bin/make install; \
    cd ${SRC_FLDR}; \
    /usr/bin/rm -rf ost; \
    # make the compounds.chemlib interchangeable: move to its own dir and link
    /usr/bin/mkdir /compound_lib; \
    /usr/bin/mv /usr/local/share/openstructure/compounds.chemlib \
                /compound_lib/; \
    ln -s /compound_lib/compounds.chemlib /usr/local/share/openstructure/; \
    # \
    # Cleanup packages only needed to compile/ install things \
    # \
    /usr/bin/yum --assumeyes remove fftw-devel gcc gcc-c++ glibc-devel \
                             libpng-devel make python3-devel \
                             sqlite-devel zlib-devel; \
    /usr/bin/yum --assumeyes autoremove; \
    /usr/bin/yum --assumeyes clean all; \
    /usr/bin/rm -rf /var/cache/yum; \
    /usr/bin/rm -rf /tmp/*

#
# Install QMEAN
#
# Version can be switched via --build-arg QMEAN_VERSION="<TAG|BRANCH>"
ARG QMEAN_VERSION="4.2.0"
ENV VERSION_QMEAN=$QMEAN_VERSION
RUN set -eo pipefail; \
    # \
    # Install dependencies for building software \
    # \
    /usr/bin/yum makecache fast; \
    /usr/bin/yum -y install epel-release; \
    /usr/bin/yum update --assumeyes; \
    /usr/bin/yum install --assumeyes gcc-c++ make python3-devel; \
    /usr/bin/yum --assumeyes clean all; \
    /usr/bin/rm -rf /var/cache/yum; \
    /usr/bin/curl -o qmean.tar.gz -L ${SCCRE_GTLB}/QMEAN/repository/${QMEAN_VERSION}/archive.tar.gz; \
    /usr/bin/mkdir qmean; \
    /usr/bin/tar --strip-components=1 -C qmean -xf qmean.tar.gz; \
    /usr/bin/rm -f qmean.tar.gz; \
    /usr/bin/mkdir qmean/build; \
    cd qmean/build; \
    /usr/local/bin/cmake .. -DOPTIMIZE=1 -DOST_ROOT=/usr/local; \
    /usr/bin/make -j ${CPUS_FOR_MAKE}; \
    /usr/bin/make check; \
    /usr/bin/make install; \
    cd ${SRC_FLDR}; \
    /usr/bin/rm -rf qmean; \
    # \
    # Cleanup packages only needed to compile/ install things \
    # \
    /usr/bin/yum --assumeyes remove gcc gcc-c++ glibc-devel make python3-devel; \
    /usr/bin/yum --assumeyes autoremove; \
    /usr/bin/yum --assumeyes clean all; \
    /usr/bin/rm -rf /var/cache/yum; \
    /usr/bin/rm -rf /tmp/*; \
    # removing python3-devel also removes the entire Python 3 installation...
    /usr/bin/yum makecache fast; \
    /usr/bin/yum -y install epel-release; \
    /usr/bin/yum install --assumeyes python3; \
    /usr/bin/yum --assumeyes clean all; \
    /usr/bin/rm -rf /var/cache/yum

# only for the hackathon
#RUN set -eo pipefail; \
#    /usr/bin/python3 -m pip install --upgrade --no-cache-dir pip; \
#    /usr/bin/python3 -m pip install --no-cache-dir pylint black


# Labels
LABEL license=""
LABEL org.openstructure.ost.license="LGPL Version 3"
LABEL org.openstructure.qmean.license="Apache License Version 2.0, January 2004"
LABEL org.openstructure.version="${OST_VERSION}"
LABEL org.openstructure.qmean.version="${QMEAN_VERSION}"
LABEL org.openstructure.ost.version.is-production="yes"
LABEL org.openstructure.qmean.version.is-production="yes"
#   This is a counter to be increased with every built & uploaded image
LABEL org.openstructure.qmean.docker.image.increment="2"
LABEL maintainer="Stefan Bienert <stefan.bienert@unibas.ch>"
LABEL vendor1="Schwede Group (schwedelab.org)"
LABEL vendor2="SIB - Swiss Institute of Bioinformatics"
LABEL vendor3="Biozentrum - University of Basel"

# create user & group to run QMEAN
# IDs can be switched via --build-arg QMEAN_USER_ID=<USER ID>,
# e.g. --build-arg QMEAN_USER_ID=$(id -u) for your own user. That would mean
# that files and directories outside of the container would belong to you and
# directories you create to be mounted as volumes inside the container are owned
# by the QMEAN user. Just make sure to match the user ID with the one running the
# container in production or make corresponding directories writable. Treat the
# group ID similar.
ARG QMEAN_USER_ID=1000
ARG QMEAN_GROUP_ID=1000
ENV QMEAN_USER="qmean" \
    QMEAN_GROUP="qmean"
RUN groupadd --gid $QMEAN_GROUP_ID $QMEAN_GROUP; \
    useradd --no-log-init -u $QMEAN_USER_ID -r -ms /bin/bash -g $QMEAN_GROUP \
            $QMEAN_USER; \
    /usr/bin/mkdir /qmean; \
    /usr/bin/mkdir /data; \
    /usr/bin/chgrp -R $QMEAN_GROUP /data; \
    /usr/bin/chown -R $QMEAN_USER /data

# Copy script to run QMEANDisCo
COPY run-qmeandisco.py /qmean/
COPY predict_ss_sa.pl /qmean/
COPY homology_sa.pl /qmean/
COPY getACCAligns.pl /qmean/

ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1

WORKDIR /data

COPY docker-entrypoint.sh /qmean
ENTRYPOINT ["/qmean/docker-entrypoint.sh"]
CMD ["qmeandisco"]

USER $QMEAN_USER

# TODO:
# - note in README how to rebuild with a dedicated user, also test that this
#   works
# - turn the run-qmeandisco.py into a proper QMEAN DisCo
#   - requires HHBlits, ACCPro and maybe more software in this image
# - in the README note how to switch to a recent compounds.chemlib
# - to be discussed: provide weekly download of compounds.chemlib?
# - copy runscript from cloned repo, not build context!
# - run as user qmean, make sure that swapping the compounds library also works
# - put into docs how to create a compounds library with this container
# - mention in docs how to build a container with dedicated uid to fit into your
#   system
# - place docker-entrypoint.sh and qmeandisco script in a non-writable folder for
#   security reasons, folder should also be not mountable
# - From Docker best practices (not possible right now, maybe they revert that
#   advice in the future):
#   - run NTP client (https://docs.docker.com/develop/dev-best-practices/)
# - debian slim:
#   - QMEAN is a Python app, start with a debian-slim Python3 base image?
#   - for base image, consider what is persistantly installed as base

# LocalWords:  QMEAN arg
