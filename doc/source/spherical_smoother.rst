Spherical Smoother
================================================================================

.. currentmodule:: qmean

Internally a structure specific graph like visitor pattern gets 
generated at the smoother objects initialization. Every position gets 
represented by a node. These nodes get connected to their neighbouring
nodes with edges, that are directly associated with their distance 
dependent weights. It is now possible to insert residue specific values in
the nodes and evaluate the whole graph. This representation gives a
speed advantage, if several scores have to be smoothed for the same
structure. QMEAN typically uses the CA position to perform smoothing.

.. literalinclude:: example_scripts/smoother_example.py

.. class:: SphericalSmoother(positions, std)

  Sets up a smoother specific for the given **positions** all smoothing
  weights get precalculated for the given **std**

  :param positions:     The positions you want to smooth
  :param std:           Controls the smoothing behaviour

  :type positions:      :class:`list` of :class:`ost.geom.Vec3`
  :type std:            :class:`float`

  .. method:: Smooth(values)

    :param values:      The values you want to smooth. They're assumed to 
                        be at the positions as defined at the smoother 
                        initialization.

    :type param:        :class:`list` of :class:`float`

    :returns:           The smoothed values
    :rtype:             :class:`list` of :class:`float`

    :raises:            :class:`RuntimeError` if size of values is inconsistent
                        with size of positions given at initialization
