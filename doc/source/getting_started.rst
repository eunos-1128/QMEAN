.. Copyright (c) 2013-2018, SIB - Swiss Institute of Bioinformatics and
.. Biozentrum - University of Basel
.. 
.. Licensed under the Apache License, Version 2.0 (the "License");
.. you may not use this file except in compliance with the License.
.. You may obtain a copy of the License at
.. 
.. http://www.apache.org/licenses/LICENSE-2.0
.. 
.. Unless required by applicable law or agreed to in writing, software
.. distributed under the License is distributed on an "AS IS" BASIS,
.. WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
.. See the License for the specific language governing permissions and
.. limitations under the License.

Getting Started
================================================================================

.. currentmodule:: qmean


QMEAN is used as a Python module either in an interactive OpenStructure shell or 
in scripts being executed with the ost OpenStructure executable. You can run all 
examples of the documentation by calling the according scripts in 
path_to_qmean/doc/source/example_scripts by typing:

  .. code-block:: bash

    ost script_name.py

Build this Documentation
--------------------------------------------------------------------------------

All information on how to build/run QMEAN can be found in the documentation.
The documentation is based on Sphinx, so you need Python, Sphinx and Make
installed on your system.

  * change to the doc directory
  * fire: "make html"
  * the documentation is built in doc/build 

Build QMEAN
--------------------------------------------------------------------------------

QMEAN requires the OpenStructure framework with all its depencies to be 
installed. Please follow the instructions on 
`openstructure.org <http://www.openstructure.org/>`_ to compile from source. 
Once you successfully compiled OpenStructure, you need following python modules:

* numpy
* scipy
* matplotlib

To compile QMEAN, change to the QMEAN directory, make a build directory and run 
CMake in there by giving it the path to your openstructure installation 
(you might want to activate optimizations with -DOPTIMIZE):

  .. code-block:: bash

    mkdir build
    cd build
    cmake ../ -DOST_ROOT=path_to_ost_build/stage -DOPTIMIZE=1

If you're working on a Debian based system, cmake might have trouble
to find the appropriate Python library. This can manually be specified.
Following additional cmake flag does the job on Ubuntu 18.04 and
Debian 9:

  .. code-block:: bash

    -DPYTHON_LIBRARIES=/usr/lib/x86_64-linux-gnu/libpython2.7.so

build:

  .. code-block:: bash
    
    make


You need the QMEAN modules to be directly accessible from within
Python. Following command enables QMEAN in your current shell:

  .. code-block:: bash

    export PYTHONPATH=path_to_qmean_build_dir/stage/lib64/python2.7/site-packages:$PYTHONPATH
    

You can add this command to .bashrc, so it is set permanently.
If everything is setup correctly, you can test the setup in an
interactive ost session by typing: from "qmean import *"
or, even better, run the example scripts from the documentation.
