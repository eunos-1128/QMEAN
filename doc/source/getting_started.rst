Getting Started
================================================================================

.. currentmodule:: qmean


QMEAN is intended to be run as a Python module. The preferred way of using it's 
functionality is from within an OpenStructure shell, or scripts being executed
with OpenStructure. As soon as you have built the documentation and QMEAN, you
can run all examples of the documentation by calling the according scripts 
lying in doc/source/example_scripts in the following way in your shell:

  .. code-block:: bash

    ost script_name.py

Building this Documentation
--------------------------------------------------------------------------------

The documentation is based on Sphinx, so you need Python with Sphinx installed,
as well as Make. In the doc directory is a Makefile, which you can fire
with "make html". The documentation should then magically appear in doc/build.

Building QMEAN
--------------------------------------------------------------------------------

QMEAN requires the OpenStructure framework with all its depencies to be 
installed. Please follow the instructions on 
`openstructure.org <http://www.openstructure.org/>`_ to compile from source. 
Do not use the precompiled libraries to avoid inconsistencies in architecture 
and used libraries (e.g. boost). Once you successfully compiled OpenStructure,
you need following software installed, with a correctly set system path.

* `DSSP <http://swift.cmbi.ru.nl/gv/dssp/>`_ 

In case you want to use QMEANBrane, you also need:

* `NACCESS <http://www.bioinf.manchester.ac.uk/naccess/>`_
* `MSMS <http://mgltools.scripps.edu/packages/MSMS/>`_

Additionally, following Python modules are required:

* Numpy
* Scipy
* Matplotlib

To compile QMEAN, change to the QMEAN directory and run CMake
by giving it the path to your openstructure installation 
(you might also want to activate the optimize flag):

  .. code-block:: bash
    
    cmake . -DOST_ROOT=path_to_ost/stage -DOPTIMIZE=1

If your working on a Debian based system it might be possible,
that the Python libraries are not found. Before you start googling
around try following additional CMake flag:

  .. code-block:: bash

    -DPYTHON_LIBRARIES=/usr/lib/x86_64-linux-gnu/libpython2.7.so

Once CMake runs through fine you can fire Make:

  .. code-block:: bash
    
    make


You need the QMEAN modules to be directly accessible from within
Python. To enable this in your current shell type:

  .. code-block:: bash

    export PYTHONPATH=path_to_qmean/stage/lib64/qmean:$PYTHONPATH
    

Alternatively you can also add this command to .bashrc, so it is
set permanently.
If everything is setup correctly, you can test the setup in an
interactive ost session by typing: from qmean import *




