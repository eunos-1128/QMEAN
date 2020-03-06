# OpenStructure image, dedicated for SMNG.

FROM centos:7.4.1708

# Arguments, optional
ARG CPUS_FOR_MAKE=1

# environment variables
ENV SCCRE_GTLB="https://git.scicore.unibas.ch/schwede" \
    SRC_FLDR="/tmp" \
    CFLAGS="-march=x86-64 -mtune=generic -isystem /usr/local/include/boost/ -isystem /usr/local/openmm/include" \
    CXXFLAGS="-march=x86-64 -mtune=generic -isystem /usr/local/include/boost/ -isystem /usr/local/openmm/include" \
    PYTHONPATH="/usr/local/lib64/python3.6/site-packages/:${PYTHONPATH}" \
    # Versions of software used, prefixed with VERSION so $VERSION_<TAB><TAB>
    # in the shell produces a list of possible software versions ;)
    VERSION_EIGEN="3.3.1" \
    VERSION_SWIG="3.0.12" \
    VERSION_BOOST="1_68_0" \
    VERSION_OPENSTRUCTURE=$OST_VERSION \
    # Define temporary dirs for the build process
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
RUN set -eo pipefail; \
    # \
    # Install dependencies for building software \
    # \
    yum makecache fast;            \
    yum -y install epel-release;   \
    yum update --assumeyes;        \
    yum install --assumeyes bzip2 bzip2-devel cmake \
    doxygen                        \
    fftw-devel                     \
    freetype-devel                 \
    gcc                            \
    gcc-c++                        \
    gcc-gfortran                   \
    gd-devel                       \
    git                            \
    git-lfs                        \
    glibc-devel                    \
    graphviz                       \
    libicu-devel                   \
    libpng-devel                   \
    libpng-static                  \
    libstdc++-static               \
    libtiff-devel                  \
    make                           \
    mariadb-devel                  \
    pcre-devel                     \
    perl                            \
    python3                         \
    python3-devel                   \
    python3-pip                     \
    readline-devel                  \
    sqlite-devel                    \
    subversion                      \
    tcsh                            \
    wget which zlib-devel zlib-static; \
    yum --assumeyes clean all;      \
    /usr/bin/rm -rf /var/cache/yum; \
    # \
    # Install Python modules \
    # \
    /usr/bin/python3 -m pip install --upgrade --no-cache-dir pip; \
    /usr/bin/python3 -m pip install --no-cache-dir numpy==1.17.4 \
                                                   scipy==1.4.0 \
                                                   matplotlib==3.1.2; \
    # \
    # Install: Eigen \
    # \
    cd ${SRC_FLDR}; \
    /usr/bin/wget -nc -O eigen.tar.bz2 \
                http://bitbucket.org/eigen/eigen/get/${VERSION_EIGEN}.tar.bz2; \
    /usr/bin/mkdir eigen/; \
    /usr/bin/tar --strip-components=1 -C eigen/ -xf eigen.tar.bz2; \
    /usr/bin/rm -rf eigen.tar.bz2; \
    /usr/bin/mv eigen/Eigen /usr/local/include; \
    /usr/bin/rm -rf eigen; \
    # \
    # Install: SWIG \
    # \
    cd ${SRC_FLDR}; \
    /usr/bin/wget -nc -O swig.tar.gz \
    https://sourceforge.net/projects/swig/files/swig/swig-${VERSION_SWIG}/swig-${VERSION_SWIG}.tar.gz; \
    /usr/bin/mkdir swig; \
    /usr/bin/tar --strip-components=1 -C swig/ -xf swig.tar.gz; \
    /usr/bin/rm -rf swig.tar.gz; \
    cd swig/; \
    ./configure --prefix=/usr/local/; \
    /usr/bin/make -j ${CPUS_FOR_MAKE}; \
    /usr/bin/make install; \
    cd ${SRC_FLDR}; \
    /usr/bin/rm -rf swig; \
    # \
    # Install: Boost \
    # \
    cd ${SRC_FLDR} && \
    /usr/bin/wget -nc -O boost.tar.bz2 \
    https://sourceforge.net/projects/boost/files/boost/1.68.0/boost_${VERSION_BOOST}.tar.bz2/download; \
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
    # Cleanup packages only needed to compile/ install things to keep the \
    # image size low. \
    # \
    yum --assumeyes remove         \
    bzip2                          \
    bzip2-devel                    \
    cmake                          \
    doxygen                        \
    fftw-devel                     \
    freetype-devel                 \
    gcc                            \
    gcc-c++                        \
    gcc-gfortran                   \
    gd-devel                       \
    git                            \
    git-lfs                        \
    glibc-devel                    \
    graphviz                       \
    libicu-devel                   \
    libpng-devel                   \
    libpng-static                  \
    libstdc++-static               \
    make                           \
    mariadb-devel                  \
    pcre-devel                     \
    perl                           \
    python3-devel                  \
    python3-pip                    \
    readline-devel                 \
    sqlite-devel                   \
    subversion                     \
    tcsh                           \
    wget                           \
    which                          \
    zlib-devel                     \
    zlib-static; \
    yum --assumeyes clean all; \
    /usr/bin/rm -rf /var/cache/yum; \
    /usr/bin/rm -rf /tmp/*; \
    # removing pip alsp removes the entire Python 3 installation...
    yum makecache fast; \
    yum -y install epel-release; \
    yum install --assumeyes python3; \
    yum --assumeyes clean all; \
    /usr/bin/rm -rf /var/cache/yum

#
# Install OpenStructure
#
# Version can be switched via --build-arg OST_VERSION="<TAG|BRANCH>"
ARG OST_VERSION="2.0.0-rc"

RUN set -eo pipefail; \
    # \
    # Install dependencies for building software \
    # \
    yum makecache fast;            \
    yum -y install epel-release;   \
    yum update --assumeyes;        \
    yum install --assumeyes        \
    bzip2                          \
    bzip2-devel                    \
    cmake                          \
    doxygen                        \
    fftw-devel                     \
    freetype-devel                 \
    gcc                            \
    gcc-c++                        \
    gcc-gfortran                   \
    gd-devel                       \
    git                            \
    git-lfs                        \
    glibc-devel                    \
    graphviz                       \
    libicu-devel                   \
    libpng-devel                   \
    libpng-static                  \
    libstdc++-static               \
    libtiff-devel                  \
    make                           \
    mariadb-devel                  \
    pcre-devel                     \
    perl                            \
    python3                         \
    python3-devel                   \
    python3-pip                     \
    readline-devel                  \
    sqlite-devel                    \
    subversion                      \
    tcsh                            \
    wget                            \
    which                           \
    zlib-devel                      \
    zlib-static;                    \
    yum --assumeyes clean all;      \
    /usr/bin/rm -rf /var/cache/yum; \
    cd ${SRC_FLDR}; \
    /usr/bin/mkdir ost; \
    /usr/bin/wget -nc -O ost.tar.gz \
         ${SCCRE_GTLB}/openstructure/repository/${OST_VERSION}/archive.tar.gz; \
    /usr/bin/tar --strip-components=1 -C ost -xf ost.tar.gz; \
    /usr/bin/rm -f ost.tar.gz; \
    /usr/bin/mkdir ost/build; \
    cd ost/build; \
    # The following to calls to `sed` are needed since the CMake version in
    # CentOS 7 is too low.
    /usr/bin/sed -i 's/cmake_minimum_required(VERSION 3.10.2 FATAL_ERROR)/cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)/' ../CMakeLists.txt; \
    /usr/bin/sed -i 's/cmake_policy(SET CMP0060 NEW)//' ../CMakeLists.txt; \
    /usr/bin/cmake ..  -DCOMPILE_TMTOOLS=OFF \
                       -DENABLE_INFO=OFF     \
                       -DENABLE_GFX=OFF      \
                       -DENABLE_GUI=OFF      \
                       -DUSE_NUMPY=OFF         \
                       -DUSE_RPATH=1         \
                       -DOPTIMIZE=1;         \
    # 1st run: get everything to build a compounds library
    /usr/bin/make chemdict_tool; \
    # get the compounds dictionary from RCSB PDB
    /usr/bin/wget -nc -O components.cif.gz \
                  ftp://ftp.wwpdb.org/pub/pdb/data/monomers/components.cif.gz; \
    stage/bin/chemdict_tool create components.cif.gz compounds.chemlib pdb; \
    stage/bin/chemdict_tool update modules/conop/data/charmm.cif compounds.chemlib charmm; \
    ## 2nd run: build OST with a compounds library
    /usr/bin/cmake .. -DCOMPOUND_LIB=compounds.chemlib; \
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
    yum --assumeyes remove bzip2 bzip2-devel cmake doxygen fftw-devel \
                           freetype-devel gcc gcc-c++ gcc-gfortran gd-devel \
                           git git-lfs glibc-devel graphviz libicu-devel \
                           libpng-devel libpng-static libstdc++-static make \
                           mariadb-devel pcre-devel perl python3-devel \
                           python3-pip readline-devel sqlite-devel subversion \
                           tcsh wget which zlib-devel zlib-static; \
    yum --assumeyes clean all; \
    /usr/bin/rm -rf /var/cache/yum; \
    /usr/bin/rm -rf /tmp/*; \
    # removing pip alsp removes the entire Python 3 installation...
    yum makecache fast; \
    yum -y install epel-release; \
    yum install --assumeyes python3; \
    yum --assumeyes clean all; \
    /usr/bin/rm -rf /var/cache/yum

#
# Install QMEAN
#
# Version can be switched via --build-arg QMEAN_VERSION="<TAG|BRANCH>"
ARG QMEAN_VERSION="4.0.0-rc"

RUN set -eo pipefail; \
    # \
    # Install dependencies for building software \
    # \
    yum makecache fast;            \
    yum -y install epel-release;   \
    yum update --assumeyes;        \
    yum install --assumeyes bzip2 bzip2-devel cmake doxygen fftw-devel \
                            freetype-devel gcc gcc-c++ gcc-gfortran gd-devel \
                            git git-lfs glibc-devel graphviz libicu-devel \
                            libpng-devel libpng-static libstdc++-static \
                            libtiff-devel make mariadb-devel pcre-devel perl \
                            python3 python3-devel python3-pip readline-devel \
                            sqlite-devel subversion tcsh wget which zlib-devel \
                            zlib-static; \
    yum --assumeyes clean all;      \
    /usr/bin/rm -rf /var/cache/yum; \
    /usr/bin/wget -nc -O qmean.tar.gz \
            ${SCCRE_GTLB}/QMEAN/repository/${QMEAN_VERSION}/archive.tar.gz; \
    /usr/bin/mkdir qmean; \
    /usr/bin/tar --strip-components=1 -C qmean -xf qmean.tar.gz; \
    /usr/bin/rm -f qmean.tar.gz; \
    /usr/bin/mkdir qmean/build; \
    cd qmean/build; \
    # The following to calls to `sed` are needed since the CMake version in
    # CentOS 7 is too low.
    /usr/bin/sed -i 's/cmake_minimum_required(VERSION 3.12.1 FATAL_ERROR)/cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)/' ../CMakeLists.txt; \
    /usr/bin/sed -i 's/cmake_policy(SET CMP0060 NEW)//' ../CMakeLists.txt; \
    /usr/bin/cmake .. -DOPTIMIZE=1 -DOST_ROOT=/usr/local; \
    /usr/bin/make -j ${CPUS_FOR_MAKE}; \
    /usr/bin/make check; \
    /usr/bin/make install; \
    cd ${SRC_FLDR}; \
    /usr/bin/rm -rf qmean; \
    # \
    # Cleanup packages only needed to compile/ install things \
    # \
    yum --assumeyes remove bzip2 bzip2-devel cmake doxygen fftw-devel \
                           freetype-devel gcc gcc-c++ gcc-gfortran gd-devel \
                           git git-lfs glibc-devel graphviz libicu-devel \
                           libpng-devel libpng-static libstdc++-static make \
                           mariadb-devel pcre-devel perl python3-devel \
                           python3-pip readline-devel sqlite-devel subversion \
                           tcsh wget which zlib-devel zlib-static; \
    yum --assumeyes clean all; \
    /usr/bin/rm -rf /var/cache/yum; \
    /usr/bin/rm -rf /tmp/*; \
    # removing pip alsp removes the entire Python 3 installation...
    yum makecache fast; \
    yum -y install epel-release; \
    yum install --assumeyes python3; \
    yum --assumeyes clean all; \
    /usr/bin/rm -rf /var/cache/yum

# Copy script to run QMEANDisCo
# Set CMD or Entrypoint?

# Labels
LABEL license=""
LABEL org.openstructure.ost.license="LGPL Version 3"
LABEL org.openstructure.qmean.license="Apache License Version 2.0, January 2004"
LABEL org.openstructure.version="${OST_VERSION}"
LABEL org.openstructure.qmean.version="${QMEAN_VERSION}"
LABEL org.openstructure.ost.version.is-production="no"
LABEL org.openstructure.qmean.version.is-production="no"
LABEL maintainer="Stefan Bienert <stefan.bienert@unibas.ch>"
LABEL vendor1="Schwede Group (schwedelab.org)"
LABEL vendor2="SIB - Swiss Institute of Bioinformatics"
LABEL vendor3="Biozentrum - University of Basel"

RUN /usr/bin/mkdir /qmean

WORKDIR /qmean

COPY run-qmeandisco.py /qmean/
