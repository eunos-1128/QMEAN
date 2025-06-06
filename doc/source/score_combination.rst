.. Copyright (c) 2013-2020, SIB - Swiss Institute of Bioinformatics and
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

Combining Scores
================================================================================

.. currentmodule:: qmean


QMEAN tries to give accurate quality estimates on global and local per
residue scale by combining different aspects / scores.
Either this happens with a simple linear combination or by more 
sophisticated multi-layer perceptrons.



Linear Score combination
--------------------------------------------------------------------------------


.. literalinclude:: example_scripts/local_scorer_example.py

Following classes allow to train and apply linear weights for arbitrary 
linear combinations of scores. The training procedure relies on the
usage of the ost table class.

.. class:: LocalScorer()

  QMEAN can combine local scores in a linear manner. Due to amino acid specific
  biases, such linear combinations are generated for all 20 standard amino 
  acids. Another problem are single values that can be NaN. An example for that
  are the first and last two residues in a torsion potential due to invalid
  dihedral angles. To solve this second problem, linear combinations are
  calculated for all possible combinations of the features. If a score
  of a particular residue with certain NaN scores has to be estimated,
  only the subset of single scores not being NaN are considered.

  .. method:: Save(filename)

    Save down the scorer using the pickle functionality of Python

    :param filename:    The filename you want to dump the 
                        :class:`LocalScorer` to
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


.. class:: GlobalScorer()

  QMEAN can combine global scores in a linear manner, exactly the same way as
  local scores. A problem are single
  values, that can be NaN. To solve this problem, linear combinations
  are calculated for all possible combinations of the features. If a score
  of a particular structure with certain NaN scores has to be estimated,
  only the subset of single scores not being NaN are considered.

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

Multi-Layer Perceptron scoring
--------------------------------------------------------------------------------

To be independent from any machine learning library, QMEAN comes with a 
lightweight multi-layer perceptron implementation: the :class:`Regressor`. 
The idea is to use whatever machine learning library to train such a network and 
transform it into a :class:`Regressor`. For the simplest training cases 
:func:`TrainRegressor` should be sufficient.

.. class:: Regressor(n_input_features, mean=None, std=None)

  Lighweight implementation of a fully connected multi-layer perceptron with 
  in-built data normalization. Before traversing the layers, data is normalized 
  as (input-mean)/std.

  :param n_input_features: The number of input features to expect
  :param mean: Mean value for every input feature (*n_input_features* values) 
               for data normalization. All values set to 0.0 if not provided.
  :param std: Standard deviation for every input feature (*n_input_features* 
              values) for data normalization. All values set to 1.0 if not 
              provided.

  :type n_input_features: :class:`int`
  :type mean: :class:`list` of :class:`float` or :class:`numpy.ndarray`
  :type mean: :class:`list` of :class:`float` or :class:`numpy.ndarray`
   

  .. method:: AddLayer(weights, activation_function, bias=None)

    Adds another layer *i*. Given the weight matrix M, the new layer
    is derived from a simple matrix multiplication M*layer[*i*-1], where
    the input layer is considered to be layer 0. Postprocessing 
    consists of adding bias values and perform a final activation.

    :param weights:  Weight matrix with shape (m,n). n is the number of
                     elements of the previous layer and m determines the
                     number of elements in the added layer.
    :param activation_function: 0: no activation, 1: ReLU (Rectified Linear 
                                Unit) activation
    :param bias:     Bias that is applied to layer before executing the 
                     activation function. Size must be consistent with
                     number of rows (m) in *weights*. All values set to 0.0 if 
                     not provided.
    :type weights:  :class:`numpy.ndarray`
    :type activation_function: :class:`int`
    :type bias:  :class:`list` of :class:`float` or :class:`numpy.ndarray`

  .. method:: Save(filepath)

    Dumps regressor in binary format

    :param filepath: Path to dump regressor
    :type filepath: :class:`str`

  .. method:: Load(filepath)

    Static method to load previously dumped regressor

    :param filepath: Path to dumped regressor
    :type filepath: :class:`str`
    
  .. method:: Predict(input)

    Normalizes input, traverses all layers an returns first element of last 
    layer: the prediction

    :param input:  Input of size *n_input_features* to feed into network
    :type input:  :class:`list` of :class:`float` or :class:`numpy.ndarray`


.. method:: TrainRegressor(data_frame, features, target, loss_function,\
                           optimizer, topology, epochs, batch_size,\
                           randomize=False)

  Trains a :class:`Regressor` using keras given a pandas dataframe as input.
  This obviously adds pandas and keras as dependency. Check the example script!

  :param data_frame:  Data with a column for each feature as well as the target
  :param features:    Describes input features, features must be present in 
                      *data_frame*
  :param target:      The target, must be present in *data_frame*
  :param loss_function: Any function that keras understands. E.g. 
                      "mean_squared_error" or "mean_absolute_error"
  :param optimizer:   Any optimizer that keras understands. E.g.
                      "sgd", "rmsprop", "adagrad", "adadelta", "adam"
  :param topology:    Every element describes number of nodes of a layer
                      with first layer being the input layer that must have 
                      the same size as *features*. The last layer is the
                      output layer that must be of size 1.
  :param epochs:      Number of training epochs
  :param batch_size:  Training batch size
  :param randomize:   Whether to randomize the order of training data prior to
                      training

  :type data_frame:   :class:`pandas.DataFrame`
  :type features:     :class:`list` of :class:`str`
  :type target:       :class:`str`
  :type loss_function: :class:`str`
  :type optimizer:    :class:`str`
  :type topology:     :class:`list` of :class:`int`
  :type epochs:       :class:`int`
  :type batch_size:   :class:`int`
  :type randomize:    :class:`bool`
  
.. literalinclude:: example_scripts/regressor_training.py


The :class:`Regressor` requires a fixed set of input scores and cannot 
flexibly adapt if certain features are missing. Examples include scores for the
first/last two residues in the :class:`qmean.TorsionPotential` due to missing
dihedral angles. The :class:`NNScorer` offers a naive solution and selects the
appropriate :class:`Regressor` in a set of alternatives that cover the possible
combinations of input scores.

.. class:: NNScorer(path)

  Organizes multiple :class:`Regressor` which cover alternative combinations of 
  input scores and selects the right one for scoring. The object must be 
  constructed manually and is loaded from disk.

  :param path: Directory containing a file 'feature_groups.json' and several 
               stored :class:`Regressor` with naming 'nn_<idx>.dat'. The first 
               is a json file containing a list of items. Each item at location 
               idx is a list of score names that represent the input for the
               the :class:`Regressor` named 'nn_<idx>.dat'.
  :type path:  :class:`str`


  .. method:: GetScore(score_dict, olc=None)

    Extracts all valid scores from *score_dict*, selects the right internal 
    :class:`NNScorer` and returns a score.

    :param score_dict: Keys relate to the internal feature group list. All 
                       values !(None || NaN) are selected and the 
                       internal list of feature groups is iterated until 
                       a group can be identified for which all values are 
                       valid. Ordering of the internal feature groups therefore 
                       matters. The values are ordered as defined in the found 
                       feature group and passed to the according 
                       :class:`Regressor`. The function returns 0.0 if no 
                       feature group can be identified.
    :param olc:        A common use case for the :class:`NNScorer` is to 
                       estimate residue specific scores. When defining *olc*,
                       a one-hot array with 20 elements is added in front
                       of the values that are passed to the selected
                       :class:`Regressor`. All elements are zero, the location
                       of *olc* in the string "ACDEFGHIKLMNPQRSTVWY" is set
                       to one. The function returns 0.0 if *olc* is not found in
                       the specified string.
    :type score_dict:  :class:`dict`
    :type olc:         :class:`str`





  
      




