Container Classes
================================================================================

.. currentmodule:: qmean

Dealing with many different statistic or potential objects can be tedious and
often result in many files on your harddrive. QMEAN offers you container
classes to gather statistic or potential objects. They can then be 
saved/loaded from/to disk at once and easily be accessed in a 
dictionary like manner.


.. class:: StatisticContainer

  The statistic container holds as many statistic objects as you want and
  allows to save them to disk in one single file.

  .. method:: Load(filename)
    
    Loads and returns a statistic container from disk

    :param filename:    Name of file
    :type filename:     :class:`str`

    :returns:           The loaded container
    :rtype:             :class:`StatisticContainer`

    :raises:            :class:`RuntimeError` if file does not exist

  .. method:: Save(filename)

    Saves a statistic container down to disk

    :param filename:    Name of file
    :type filename:     :class:`str`


  .. method:: GetKeys()

    :returns:           All keys of the statistics present in the container
    :rtype:             :class:`list` of :class:`str`

  .. method:: __getitem__(name)

    :param name:        Name of the requested statistic
    :type name:         :class:`str`
    :returns:           The requested statistic object
    :rtype:             Specific instance of :class:`StatisticBase`
    :raises:            :class:`RuntimeError` if statistic object is not present

  .. method:: __setitem__(name, statistic)

    :param statistic:   Statistic to be added
    :param name:        Name/key of the statistic to be added
    :type name:         :class:`str`
    :type statistic:    :class:`StatisticBase`

  .. method:: __delitem__(name)

    :param name:        Statistic to be removed
    :type name:         :class:`str`

  .. method:: __contains__(name)

    :param name:        Name of statistic to be checked for presence
    :type name:         :class:`str`
    :returns:           Whether requested statistic is present in container
    :rtype:             :class:`bool`

  .. method:: __len__()

    :returns:           Number of statistics in container
    :rtype:             :class:`int`





.. class:: PotentialContainer

  The potential container holds as many potential objects as you want and
  allows to save them to disk in one single file.

  .. method:: Load(filename)
    
    Loads and returns a potential container from disk

    :param filename:    Name of file
    :type filename:     :class:`str`

    :returns:           The loaded container
    :rtype:             :class:`PotentialContainer`

    :raises:            :class:`RuntimeError` if file does not exist

  .. method:: Save(filename)

    Saves a potential container down to disk

    :param filename:    Name of file
    :type filename:     :class:`str`


  .. method:: GetKeys()

    :returns:           All keys of the potentials present in the container
    :rtype:             :class:`list` of :class:`str`

  .. method:: __getitem__(name)

    :param name:        Name of the requested potential
    :type name:         :class:`str`
    :returns:           The requested potential object
    :rtype:             Specific instance of :class:`PotentialBase`
    :raises:            :class:`RuntimeError` if potential object is not present

  .. method:: __setitem__(name, potential)

    :param potential:   Potential to be added
    :param name:        Name/key of the potential to be added
    :type name:         :class:`str`
    :type potential:    :class:`PotentialBase`

  .. method:: __delitem__(name)

    :param name:        Potential to be removed
    :type name:         :class:`str`

  .. method:: __contains__(name)

    :param name:        Name of potential to be checked for presence
    :type name:         :class:`str`
    :returns:           Whether requested potential is present in container
    :rtype:             :class:`bool`

  .. method:: __len__()

    :returns:           Number of potentials in container
    :rtype:             :class:`int`

