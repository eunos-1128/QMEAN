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

Membrane Location
================================================================================

An implicit solvation model can be used to detect the transmembrane region of
a membrane protein. QMEANBrane uses an implementation described by the 
Orientation of Proteins in Membrane (OPM) database [lomize2006]_.
The actual method is implemented in :meth:`ost.mol.alg.FindMembrane`.


.. [lomize2006] Lomize AL, Pogozheva ID, Lomize MA, Mosberg HI (2006). Positioning of proteins in membranes: A computational approach. Protein Science 15, 1318-1333

