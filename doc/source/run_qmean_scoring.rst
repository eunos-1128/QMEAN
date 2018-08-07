Run the QMEAN Scoring Functions
================================================================================


.. currentmodule:: qmean


There are two conventient functions applying the quality assessment to
protein structure models. One general function and another one specifically
tailered towards the quality assessment of membrane protein models as
described in QMEANBrane.



Assessing the quality for soluble protein models
--------------------------------------------------------------------------------


.. literalinclude:: example_scripts/assess_model_quality_example.py

.. method:: AssessModelQuality(model,[output_dir='.', plots=True,  \
                               local_scores=True, global_scores=True \
                               table_format="ost",psipred=None, dc=None \
                               accpro=None, assign_bfactors=True, \
                               settings=None])

  :param model:         The model you want to assess
  :param output_dir:    The directory where to write the output 
  :param plots:         Whether you want to create plots of the output
  :param local_scores:  Whether you want to generate local scores
  :param global_scores: Whether you want to generate global scores
  :param table_format:  The format you want to save the results in, can
                        be one of

                        * "ost"     ost-specific format (human_readable)
                        * "csv"     comma separated values (human readable)
                        * "pickle"  pickled byte stream (binary)
                        * "html"    HTML table

  :param psipred:       Information from the PSIPRED prediction, QMEAN will run
                        through without it being set but for optimal performance
                        this variable must be set.
  :param accpro:        Information from the ACCPRO prediction, QMEAN will run
                        through without it being set but for optimal performance
                        this variable must be set.

  :param dc:            Thats what let QMEAN evolve to QMEANDisCo. 
                                            

  :param assign_bfactors: If set to True, the local scores get assigned to the
                          **model** as bfactors

  :param settings:      Object that contains paths to potentials etc.
                        If set to None, an object with default paths is
                        generated with: settings=conf.SwissmodelSettings().
                        Must have following members:
 
                          * local_potentials: path to :class:`PotentialContainer` 
                                              applied for local scoring.
                          * global_potentials: path to :class:`PotentialContainer`
                                               applied for global scoring
                          * local_scorer: path to :class:`LocalScorer` to 
                                          linearly combine per-residue scores
                          * global_scorer: path to :class:`GlobalScorer` to 
                                           linearly combine global scores
                          * reference_tab: path to :class:`ost.Table` containing
                                           reference data.
                          * disco_tree: path to random forest, that predicts
                                        local QMEANDisCo scores. The tree is pickled
                                        and dependent on the version of sklearn
                                        it has been generated with. It's highly 
                                        recommended to retrain it on your system using
                                        the data provided in 
                                        <path_to_qmean>/data/qmean/scorer/tree_generation


  :type model:          :class:`ost.mol.EntityHandle` / :class:`ost.mol.EntityView`
  :type output_dir:     :class:`str`
  :type plots:          :class:`bool`
  :type local_scores:   :class:`bool`
  :type global_scores:  :class:`bool`
  :type table_format:   :class:`str`
  :type psipred:        :class:`PSIPREDHandler` or a :class:`list` 
                        thereof (one element per chain)
  :type accpro:         :class:`ACCPROHandler` or a :class:`list`
                        thereof (one element per chain)
  :type dc:             :class:`DisCoContainer` or a :class:`list`
                        thereof
  :type assign_bfactors: :class:`bool`



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
                        When set to None, the FindMembrane function will
                        be called with **model** as input.

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
  :type mem_param:      :class:`FindMemParam`
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


