Membrane Location
================================================================================

An implicit solvation model can be used to detect the transmembrane region of
a membrane protein. QMEANBrane uses an implementation described by the 
Orientation of Proteins in Membrane (OPM) database [lomize2006]_.
The **membrane_axis** gets constructed using following sequence of rotations:

* Start with an initial **axis**
* Rotate it around a **tilt_axis** by angle **tilt**
* Rotate it around **axis** by angle **angle**

Having this axis, the membrane center can be constructed by moving **pos**
Angstrom from the origin in direction of the **membrane_axis**. The predicted
transmembrane planes can then be described by using the **membrane_axis** as
normal and the positions relative to the membrane center +-1/2 **width** 
Angstrom along the **membrane_axis**.


.. literalinclude:: example_scripts/membrane_example.py


.. class:: FindMemParam

  The result object when searching for a membrane.

  .. attribute:: axis

    The initial search axis to which everything is relative to

  .. attribute:: tilt_axis

    The axis around which the membrane axis gets tilted in a first step

  .. attribute:: tilt

    The according angle in radians

  .. attribute:: angle

    The angle in radians to rotate the tilted membrane axis around the initial 
    search axis

  .. attribute:: membrane_axis

    The transformed **axis**, that describes the membrane orientation

  .. attribute:: width

    The membrane width in Angstrom

  .. attribute:: pos

    The membrane center as described in Angstrom from the origin in direction
    membrane axis

  .. attribute:: energy

    The pseudo energy returned from the implicit solvation model

  .. method:: Save(filename)

    Dumps the object down to disk

    :param filename:    The filename you want to dump the object into
    :type filename:     :class:`str`

  .. method:: Load(filename)

    Static load function to load a FindMemParam object from disk

    :param filename:    The filename you want to load
    :type filename:     :class:`str`

    :returns:           The loaded object
    :rtype:             :class:`FindMemParam`

    :raises:            :class:`RuntimeError` when file doesn't exist



.. method:: FindMembrane(ent, surf, asa)

  Detects a transmembrane region based on an implicit solvation model

  :param ent:           The entity handle from which you want to calculate
                        the transmembrane region. All hydrogens must be
                        stripped off!

  :param surf:          Surface as calculated by MSMS 

  :param asa:           Solvent accessibilities as calculated by NACCESS
                        for every atom in **ent**

  :type ent:            :class:`ost.mol.EntityHandle`
  :type surf:           :class:`ost.mol.SurfaceHandle`
  :type asa:            :class:`list` of :class:`float`

  :returns:             Parameters, that define the predicted transmembrane region
  :rtype:               :class:`FindMemParam`

  :raises:              :class:`RuntimeError` if number of atoms in **ent** and size
                        of **asa** is inconsistent.



.. [lomize2006] Lomize AL, Pogozheva ID, Lomize MA, Mosberg HI (2006). Positioning of proteins in membranes: A computational approach. Protein Science 15, 1318-1333




