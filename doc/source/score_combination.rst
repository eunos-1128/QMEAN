Combining Scores
================================================================================

.. currentmodule:: qmean


Following classes allow to train and apply linear weights for arbitrary 
linear combinations of scores. The training procedure relies on the
usage of the ost table class.




Combining Single Residue Scores
--------------------------------------------------------------------------------

QMEAN combines local scores in a linear manner. Due to amino acid specific
biases, such linear combinations are generated for all 20 standard amino 
acids. Another problem are single values, that can be NaN. An example for that
are the first and last two residues in a torsion potential due to invalid
dihedral angles. To solve this second problem, linear combinations get
calculated for all possible combinations of the features. If a score
of a particular residue with certain NaN scores has to be estimated,
only the subset of single scores not being NaN are considered.

.. literalinclude:: example_scripts/local_scorer_example.py

.. class:: LocalScorer()

  .. method:: Save(filename)

    Save down the scorer using the pickle functionality of Python

    :param filename:    The filename you want to dump the 
                        :class:`LocalScorer` into
    :type filename:     :class:`str`

  .. method:: Load(filename)

    Static loading function for a dumped :class:`LocalScorer`

    :param filename:    Name of file you want to load
    :type filename:     :class:`str`
    :returns:           The loaded scorer 
    :rtype:             :class:`LocalScorer`
    :raises:            :class:`RuntimeError` if filename does no exist
                        or when the thing can't be loaded   

  .. method:: TrainNewType(features, tab, residue_type, [method="leastsqr"])

    Trains Linear combinations for all possible amino acids and all possible
    combinations of the **features** based on the data given in **tab**.

    :param features:    Name of the features as they are stored in **tab**
    :param tab:         Table containing the training information. Several
                        columns are mandatory: 
                          * "AA": one letter code for data points
                          * "target": target value for data points
                          * Column for every element in **features**
    :param residue_type: Every scorer can train several score combinations
                         you might want to distinct between secondary
                         structure of whatever
    :param method:       Method to generate the linear combination.
                         possible values are:
                           * "leastsqr": uses lstsq functionality from the
                             scipy.linalg module
                           * "robust": uses the robustfit function 
                             from matlab. You need to have matlab installed
                             and a running mlabwrap (google is your friend)

    :type features:     :class:`list` of :class:`str`
    :type tab:          :class:`ost.Table`
    :type residue_type: :class:`str`
    :type method:       :class:`str`

    :raises:            :class:`RuntimeError` if one of the mandatory fields
                        is not present in **tab**, the method is invalid,
                        or you choose "robust" as method and mlabwrap is
                        not setup properly


  .. method:: GetLocalScore(residue_type, AA, data)

    :param residue_type: Type of score combination, thats the one you called
                         the TrainNewType function with
    :param AA:           The identity of the amino acid in form of a
                         one letter code
    :param data:         Dictionary with the features as keys and the
                         according data as values 
    :type residue_type: :class:`str`
    :type AA:           :class:`str`
    :type data:         :class:`dict`

    :returns:           The Linear combination for the requested features
    :rtype:             :class:`float`

    :raises:            :class:`RuntimeError` If the scorer has not been
                        trained for **residue_type** or **AA** is unknown





Combining Global Scores
--------------------------------------------------------------------------------
QMEAN combines global scores in a linear manner, exactly the same way as
local scores. A problem are single
values, that can be NaN. To solve this problem, linear combinations
get calculated for all possible combinations of the features. If a score
of a particular structure with certain NaN scores has to be estimated,
only the subset of single scores not being NaN are considered.


.. class:: GlobalScorer()

  .. method:: Save(filename)

    Save down the scorer using the pickle functionality of Python

    :param filename:    The filename you want to dump the 
                        :class:`GlobalScorer` into
    :type filename:     :class:`str`

  .. method:: Load(filename)

    Static loading function for a dumped :class:`GlobalScorer`

    :param filename:    Name of file you want to load
    :type filename:     :class:`str`
    :returns:           The loaded scorer 
    :rtype:             :class:`GlobalScorer`
    :raises:            :class:`RuntimeError` if filename does no exist
                        or when the thing can't be loaded   

  .. method:: TrainNewType(features, tab, structure_type, [method="leastsqr"])

    Trains Linear combinations for all possible amino acids and all possible
    combinations of the **features** based on the data given in **tab**.

    :param features:    Name of the features as they are stored in **tab**
    :param tab:         Table containing the training information. Several
                        columns are mandatory: 
                          * "target": target value for data points
                          * Column for every element in **features**
    :param structure_type: Every scorer can train several score combinations,
                         this string is the according identifier
    :param method:       Method to generate the linear combination.
                         possible values are:
                           * "leastsqr": uses lstsq functionality from the
                             scipy.linalg module
                           * "robust": uses the robustfit function 
                             from matlab. You need to have matlab installed
                             and a running mlabwrap (google is your friend)

    :type features:     :class:`list` of :class:`str`
    :type tab:          :class:`ost.Table`
    :type structure_type: :class:`str`
    :type method:       :class:`str`

    :raises:            :class:`RuntimeError` if one of the mandatory fields
                        is not present in **tab**, the method is invalid,
                        or you choose "robust" as method and mlabwrap is
                        not setup properly


  .. method:: GetGlobalScore(structure_type,  data)

    :param residue_type: Type of score combination, thats the one you called
                         the TrainNewType function with
    :param data:         Dictionary with the features as keys and the
                         according data as values 
    :type residue_type: :class:`str`
    :type data:         :class:`dict`

    :returns:           The Linear combination for the requested features
    :rtype:             :class:`float`

    :raises:            :class:`RuntimeError` If the scorer has not been
                        trained for **structure_type**