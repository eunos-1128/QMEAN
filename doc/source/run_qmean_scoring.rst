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

Run the QMEAN Scoring Functions
================================================================================


.. currentmodule:: qmean


We distinguish between classical QMEAN to assess quality of soluble protein 
models and the membrane specific version QMEANBrane. 

Assessing the quality for soluble protein models
--------------------------------------------------------------------------------


.. literalinclude:: example_scripts/assess_model_quality_example.py

.. class:: QMEANScorer(model,[psipred=None, accpro=None, dc=None, \
                              settings=conf.SwissmodelSettings(), \
                              use_nn = True])

  The QMEANScorer allows to calculate various scores from the QMEAN 
  universe. They can be accessed through the class attributes that 
  follow a lazy evaluation scheme, i.e. scores are only calculated
  if they're actually requested. 
  Plotting functions are also available to visualize the calculated scores.

  :param model:         The model you want to assess. Only the selection
                        returned when calling model.Select("peptide=true")
                        will be used for evaluation. Residue numbers are used
                        to refer to local per residue scores. To avoid 
                        confusion, we only allow the numeric component to be 
                        present. An error will be thrown if there's an insertion
                        code.
  :param psipred:       Information from the PSIPRED prediction, QMEAN will run
                        through without it being set but for optimal performance
                        this variable must be set.
  :param accpro:        Information from the ACCPRO prediction, QMEAN will run
                        through without it being set but for optimal performance
                        this variable must be set.
  :param dc:            Thats what let QMEAN evolve to QMEANDisCo, QMEAN will
                        run through without it being set but for optimal
                        performance this variable must be set.

  :param settings:      Contains paths to potentials etc.
                        If set to None, an object with default paths is
                        generated with: settings=conf.SwissmodelSettings().
                        Must have following members:
 
                          * local_potentials: path to :class:`PotentialContainer` 
                                              applied for local scoring.
                          * global_potentials: path to :class:`PotentialContainer`
                                               applied for global scoring 
                                               (qmean4 / qmean6)
                          * local_scorer: path to :class:`NNScorer` to estimate
                                          local scores. this is used if *use_nn*
                                          is set to True (default)
                          * linear_local_scorer: path to :class:`LocalScorer` to 
                                                 linearly combine per-residue scores.
                                                 this is only used if *use_nn* is set
                                                 to False. 
                          * global_scorer: path to :class:`GlobalScorer` to 
                                           linearly combine global scores 
                                           (qmean4 / qmean6)
                          * reference_tab: path to :class:`ost.Table` containing
                                           reference data. To calculate the global
                                           Z-scores 
  :param use_nn:        Whether to use neural network to calculate local scores. 
                        This is the prefered way of doing things. Set this to False
                        if you want to fallback to the good old times when linear
                        combinations were sufficient.


  :type model:          :class:`ost.mol.EntityHandle` / :class:`ost.mol.EntityView`
  :type psipred:        :class:`PSIPREDHandler` or a :class:`list` 
                        thereof (one element per chain)
  :type accpro:         :class:`ACCPROHandler` or a :class:`list`
                        thereof (one element per chain)
  :type dc:             :class:`DisCoContainer` or a :class:`list` thereof
  :type settings:       Custom object with specified member
  :type use_nn:         :class:`bool`


  .. method:: QMEAN4SliderPlot(out_path)

    Dumps the slider plot representing the QMEAN4 Z-Score and the Z-Scores of
    its components.

    :param out_path:    Path to dump the plot, must have image format file
                        ending, e.g. awesome_plot.png
    :type out_path:     :class:`str`


  .. method:: QMEAN6SliderPlot(out_path)

    Dumps the slider plot representing the QMEAN6 Z-Score and the Z-Scores of
    its components.

    :param out_path:    Path to dump the plot, must have image format file
                        ending, e.g. awesome_plot.png
    :type out_path:     :class:`str`


  .. method:: QMEAN4ReferencePlot(out_path)

    Dumps the reference plot with the QMEAN4 score of the *model* compared to
    the QMEAN4 scores of a large set of non-redundant X-ray structures.

    :param out_path:    Path to dump the plot, must have image format file
                        ending, e.g. awesome_plot.png
    :type out_path:     :class:`str`


  .. method:: QMEAN6ReferencePlot(out_path)

    Dumps the reference plot with the QMEAN6 score of the *model* compared to
    the QMEAN6 scores of a large set of non-redundant X-ray structures.

    :param out_path:    Path to dump the plot, must have image format file
                        ending, e.g. awesome_plot.png
    :type out_path:     :class:`str`


  .. method:: LocalProfilePlot(out_path, [chain=None])

    Dumps a line plot representing local per residue scores.

    :param out_path:    Path to dump the plot, must have image format file
                        ending, e.g. awesome_plot.png
    :param chain:       If None, the local score profiles of all chains end up
                        in the same plot. You can specify a chain name to 
                        only plot the score profile of one chain.
    :type out_path:     :class:`str`
    :type chain:        :class:`str`


  .. method:: AssignModelBFactors()

    Assigns local scores to the bfactors of the internally scored model 
    (the result of *model*.Select("peptide=true")). This operation happens in 
    place, so *model* is directly modified.


  .. method:: model

    The internally scored model (the result of 
    *model*.Select("peptide=true"))

    :rtype:             :class:`ost.mol.EntityView` 


  .. attribute:: qmean4_score

    The QMEAN4 score of *model*, i.e. four statistical potentials linearly combined 
    to get a global score.

  .. attribute:: qmean4_z_score

    The QMEAN4 Z-score of *model*, i.e. its QMEAN4 score compared to what one would 
    expect from high resolution X-ray structures

  .. attribute:: qmean4_components

    Dictionary containing scores that contribute to QMEAN4 in their normalized 
    and Z-score version. Keys of the dictionary: ["interaction", "cbeta", 
    "packing", "torsion", "interaction_z_score", "cbeta_z_score", 
    "packing_z_score", "torsion_z_score"]

  .. attribute:: qmean6_score

    The QMEAN6 score of *model*, i.e. four statistical potentials and two agreement
    terms linearly combined to get a global score.

  .. attribute:: qmean6_z_score

    The QMEAN6 Z-score of *model*, i.e. its QMEAN6 score compared to what one would 
    expect from high resolution X-ray structures

  .. attribute:: qmean6_components

    Dictionary containing scores that contribute to QMEAN6 in their normalized 
    and Z-score version. Keys of the dictionary: ["interaction", "cbeta", 
    "packing", "torsion", "ss_agreement", "acc_agreement", 
    "interaction_z_score", "cbeta_z_score", "packing_z_score", "torsion_z_score", 
    "ss_agreement_z_score", "acc_agreement_z_score"]

  .. attribute:: local_scores

    Predicted local per residue scores of *model*. This attribute returns a 
    dictionary with keys being the chain names in *model* and values another
    set of dictionaries. The key in those per chain dictionaries are the 
    residue numbers as integers (please note, that insertion codes throw an
    error if present in *model*) and values the predicted local score.

  .. attribute:: avg_local_score

    The average of the the local scores










Assessing the quality for membrane protein models
--------------------------------------------------------------------------------

.. literalinclude:: example_scripts/assess_membrane_model_quality.py

.. method:: AssessMembraneModelQuality(model, [mem_param = None, \ 
                                       output_dir='.', \
                                       plots = True, table_format='ost', \
                                       psipred=None, accpro=None, dc=None\
                                       assign_bfactors=True, settings=None])

  :param model:         The model you want to assess

  :param mem_param:     Parameters defining the membrane planes, from which
                        the core/interface membrane residues can be identified.
                        When set to None, the ost.mol.alg.FindMembrane function 
                        will be called with **model** as input.

  :param output_dir:    The directory where to write the output

  :param plots:         Whether you want to create plots of the output

  :param table_format:  The format you want to save the results in, can
                        be one of

                        * "ost"     ost-specific format (human_readable)
                        * "csv"     comma separated values (human readable)
                        * "pickle"  pickled byte stream (binary)
                        * "html"    HTML table

  :param psipred:       Information from the PSIPRED prediction, QMEAN will run
                        through without it being set but for optimal performance
                        this variable must be set. Please note, that the PSIPRED
                        prediction won't be used for assessing the quality in
                        the transmembrane part, only in the soluble part.   

  :param accpro:        Information from the ACCPRO prediction, QMEAN will run
                        through without it being set but for optimal performance
                        this variable must be set. Please note, that the ACCPRO
                        prediction won't be used for assessing the quality in
                        the transmembrane part, only in the soluble part.

  :param dc:            Information from homologous templates with known 
                        structure, QMEAN will run through without it being set but for optimal performance this variable must be set. 

  :param assign_bfactors: If set to True, the local scores get assigned to the
                          **model** as bfactors

  :param settings:      Object that contains paths to potentials etc.
                        If set to None, an object with default paths is
                        generated with: settings=conf.MembraneSettings().
                        Must have following members:
 
                          * local_potentials_soluble: path to :class:`PotentialContainer` 
                                              applied for soluble local scoring.
                          * global_potentials_soluble: path to :class:`PotentialContainer`
                                               applied for soluble global scoring
                          * local_potentials_membrane: path to :class:`PotentialContainer` 
                                              applied for local scoring in membrane and interface.
                          * global_potentials_membrane: path to :class:`PotentialContainer`
                                               applied for global scoring in membrane and interface.
                          * local_scorer_soluble: path to :class:`LocalScorer` to 
                                          linearly combine soluble per-residue scores
                          * local_scorer_membrane: path to :class:`LocalScorer` to 
                                          linearly combine membrane and interface per-residue scores
                          * global_scorer: path to :class:`GlobalScorer` to 
                                           linearly combine soluble global scores
                          * reference_tab: path to :class:`ost.Table` containing
                                           soluble reference data.


  :type model:          :class:`ost.mol.EntityHandle` / :class:`ost.mol.EntityView`
  :type mem_param:      :class:`ost.mol.alg.FindMemParam`
  :type output_dir:     :class:`str`
  :type plots:          :class:`bool`
  :type table_format:   :class:`str`
  :type psipred:        :class:`PSIPREDHandler`
  :type accpro:         :class:`ACCPROHandler`
  :type assign_bfactors: :class:`bool`


The following code is not necessarily related to the usage of 
QMEAN/QMEANBrane but rather reproduces a plot based on data produced
in the last script (Fig4 in QMEANBrane publication).

.. literalinclude:: example_scripts/reproduce_fig4_from_publication.py


