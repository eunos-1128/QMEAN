FROM centos:7.4.1708

# Arguments, optional
ARG CPUS_FOR_MAKE=1

# environment variables
ENV SCCRE_GTLB="https://git.scicore.unibas.ch/schwede" \
    SRC_FLDR="/tmp" \
    CFLAGS="-march=x86-64 -mtune=generic -isystem /usr/local/include/boost/ -isystem /usr/local/openmm/include" \
    CXXFLAGS="-march=x86-64 -mtune=generic -isystem /usr/local/include/boost/ -isystem /usr/local/openmm/include" \
    PYTHONPATH="/usr/local/lib64/python3.6/site-packages/:${PYTHONPATH}" \
    # Define temporary dirs for the build process \
    TEMP="/tmp" \
    TEMPDIR="/tmp" \
    TMPDIR="/tmp" \
    TMP="/tmp" \
    LD_LIBRARY_PATH="/usr/local/lib64:/usr/local/lib:$LD_LIBRARY_PATH"

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
ENV VERSION_EIGEN="3.3.1" \
    VERSION_SWIG="3.0.12" \
    VERSION_BOOST="1.68.0" \
    VERSION_CMAKE="3.17.3"
RUN set -eo pipefail; \
    # \
    # Install dependencies for building software \
    # \
    yum makecache fast; \
    yum -y install epel-release; \
    yum update --assumeyes; \
    yum install --assumeyes openssl-devel bzip2 gcc-c++ make pcre-devel python3 \
                            python3-devel which zlib-static; \
    yum --assumeyes clean all;      \
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
    # Install a recent version of CMake \
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
     cd ${SRC_FLDR} && \
     /usr/bin/curl -o boost.tar.bz2 -L \
     https://sourceforge.net/projects/boost/files/boost/${VERSION_BOOST}/boost_${VERSION_BOOST//\./_}.tar.bz2/download; \
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
     # Install: SWIG \
     # \
     cd ${SRC_FLDR}; \
     /usr/bin/curl -o swig.tar.gz -L \
     https://sourceforge.net/projects/swig/files/swig/swig-${VERSION_SWIG}/swig-${VERSION_SWIG}.tar.gz; \
     /usr/bin/mkdir swig; \
     /usr/bin/tar --strip-components=1 -C swig/ -xf swig.tar.gz; \
     /usr/bin/rm -rf swig.tar.gz; \
     cd swig/; \
     ./configure --without-python --prefix=/usr/local/; \
     /usr/bin/make -j ${CPUS_FOR_MAKE}; \
     /usr/bin/make install; \
     cd ${SRC_FLDR}; \
     /usr/bin/rm -rf swig; \
     # \
     # Cleanup packages only needed to compile/ install things to keep the \
     # image size low. \
     # \
     yum --assumeyes remove openssl-devel bzip2 gcc gcc-c++ glibc-devel make pcre-devel perl \
                            python3-devel python3-pip which zlib-devel \
                            zlib-static; \
     yum --assumeyes autoremove; \
     yum --assumeyes clean all; \
     /usr/bin/rm -rf /var/cache/yum; \
     /usr/bin/rm -rf /tmp/*; \
     # removing python3-devel alsp removes the entire Python 3 installation...
     yum makecache fast; \
     yum -y install epel-release; \
     yum install --assumeyes python3; \
     yum --assumeyes clean all; \
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
    yum makecache fast;            \
    yum -y install epel-release;   \
    yum update --assumeyes;        \
    yum install --assumeyes fftw fftw-devel gcc-c++ libpng \
                            libpng-devel libtiff-devel make python3-devel \
                            sqlite-devel \
    yum --assumeyes clean all;      \
    /usr/bin/rm -rf /var/cache/yum; \
    cd ${SRC_FLDR}; \
    /usr/bin/mkdir ost; \
    /usr/bin/curl -o ost.tar.gz -L \
         ${SCCRE_GTLB}/openstructure/repository/${OST_VERSION}/archive.tar.gz; \
    /usr/bin/tar --strip-components=1 -C ost -xf ost.tar.gz; \
    /usr/bin/rm -f ost.tar.gz; \
    /usr/bin/mkdir ost/build; \
    cd ost/build; \
    # The following to calls to `sed` are needed since the CMake version in
    # CentOS 7 is too low.
    #/usr/bin/sed -i 's/cmake_minimum_required(VERSION 3.12.1 FATAL_ERROR)/cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)/' ../CMakeLists.txt; \
    /usr/bin/sed -i 's/cmake_policy(SET CMP0060 NEW)//' ../CMakeLists.txt; \
    /usr/local/bin/cmake ..  -DCOMPILE_TMTOOLS=OFF \
                       -DENABLE_INFO=OFF     \
                       -DENABLE_GFX=OFF      \
                       -DENABLE_GUI=OFF      \
                       -DUSE_NUMPY=OFF         \
                       -DUSE_RPATH=1         \
                       -DOPTIMIZE=1;         \
    # 1st run: get everything to build a compounds library
    /usr/bin/make chemdict_tool; \
    # get the compounds dictionary from RCSB PDB
    /usr/bin/curl -o components.cif.gz -L \
                  ftp://ftp.wwpdb.org/pub/pdb/data/monomers/components.cif.gz; \
    stage/bin/chemdict_tool create components.cif.gz compounds.chemlib pdb; \
    stage/bin/chemdict_tool update modules/conop/data/charmm.cif compounds.chemlib charmm; \
    ## 2nd run: build OST with a compounds library
    /usr/local/bin/cmake .. -DCOMPOUND_LIB=compounds.chemlib; \
    /usr/bin/make -j ${CPUS_FOR_MAKE}; \
    /usr/bin/make check; \
    /usr/bin/make install; \
    cd ${SRC_FLDR}; \
    /usr/bin/rm -rf ost; \
    # make the compounds.chemlib interchangeable: move to its own dir and link
    /usr/bin/mkdir /compound_lib; \
    /usr/bin/mv /usr/local/share/openstructure/compounds.chemlib /compound_lib/; \
    ln -s /compound_lib/compounds.chemlib /usr/local/share/openstructure/; \
    # \
    # Cleanup packages only needed to compile/ install things \
    # \
    yum --assumeyes remove fftw-devel gcc gcc-c++ glibc-devel \
                           libpng-devel make perl python3-devel \
                           sqlite-devel zlib-devel; \
    yum --assumeyes autoremove; \
    yum --assumeyes clean all; \
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
    yum makecache fast;            \
    yum -y install epel-release;   \
    yum update --assumeyes;        \
    yum install --assumeyes gcc-c++ make python3-devel; \
    yum --assumeyes clean all;      \
    /usr/bin/rm -rf /var/cache/yum; \
    /usr/bin/curl -o qmean.tar.gz -L \
            ${SCCRE_GTLB}/QMEAN/repository/${QMEAN_VERSION}/archive.tar.gz; \
    /usr/bin/mkdir qmean; \
    /usr/bin/tar --strip-components=1 -C qmean -xf qmean.tar.gz; \
    /usr/bin/rm -f qmean.tar.gz; \
    /usr/bin/mkdir qmean/build; \
    cd qmean/build; \
    # The following to calls to `sed` are needed since the CMake version in
    # CentOS 7 is too low.
    #/usr/bin/sed -i 's/cmake_minimum_required(VERSION 3.12.1 FATAL_ERROR)/cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)/' ../CMakeLists.txt; \
    /usr/bin/sed -i 's/cmake_policy(SET CMP0060 NEW)//' ../CMakeLists.txt; \
    /usr/local/bin/cmake .. -DOPTIMIZE=1 -DOST_ROOT=/usr/local; \
    /usr/bin/make -j ${CPUS_FOR_MAKE}; \
    /usr/bin/make check; \
    /usr/bin/make install; \
    cd ${SRC_FLDR}; \
    /usr/bin/rm -rf qmean; \
    # \
    # Cleanup packages only needed to compile/ install things \
    # \
    yum --assumeyes remove gcc gcc-c++ glibc-devel make perl \
                           python3-devel; \
    yum --assumeyes autoremove; \
    yum --assumeyes clean all; \
    /usr/bin/rm -rf /var/cache/yum; \
    /usr/bin/rm -rf /tmp/* \
    # removing python3-devel also removes the entire Python 3 installation...
    yum makecache fast; \
    yum -y install epel-release; \
    yum install --assumeyes python3; \
    yum --assumeyes clean all; \
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
LABEL org.openstructure.ost.version.is-production="no"
LABEL org.openstructure.qmean.version.is-production="no"
#   This is a counter to be increased with every built & uploaded image
LABEL org.openstructure.qmean.docker.image.increment="1"
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
RUN groupadd --gid $QMEAN_GROUP_ID $QMEAN_GROUP && \
    useradd --no-log-init -u $QMEAN_USER_ID -r -ms /bin/bash -g $QMEAN_GROUP \
            $QMEAN_USER; \
    /usr/bin/mkdir /qmean; \
    /usr/bin/chgrp -R $QMEAN_GROUP /qmean; \
    /usr/bin/chown -R $QMEAN_USER /qmean

# Copy script to run QMEANDisCo
COPY run-qmeandisco.py /qmean/

ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1

WORKDIR /qmean

COPY docker-entrypoint.sh /qmean
ENTRYPOINT ["/qmean/docker-entrypoint.sh"]
CMD ["qmeandisco"]

# TODO:
# - check if SWIG is needed
# - make QMEAN DisCo the default command, make it run at PID 0 (or is it 1?
#   Check the docs!)
# - create a dedicated user (restricts writing of output because of permissions)
#   note in README how to rebuild with a dedicated user, also
#   test that this works
# - turn the run-qmeandisco.py into a proper QMEAN DisCo
#   - requires HHBlits, ACCPro and maybe more software in this image
# - in the README note how to switch to a recent compounds.chemlib
# - to be discussed: provide weekly download of compounds.chemlib?
# - copy runscript from cloned repo, not build context!

# # LocalWords:  QMEAN arg
