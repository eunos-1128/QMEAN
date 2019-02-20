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

DisCo
================================================================================

.. currentmodule:: qmean

DisCo evaluates the consistency of pairwise CA distances with an ensemble of 
constraints extracted from templates homologous to the target structure.
The following script gives some idea how the underlying object, the
:class:`DisCoContainer`, works.

.. literalinclude:: example_scripts/disco_example.py


.. class:: DisCoContainer(seqres)

  The DisCoContainer gathers constraints from an arbitrary amount of templates.
  Once all templates are added, the actual constraint functions can be 
  generated.

  :param seqres:        The SEQRES to which all templates relate
  :type seqres:         :class:`ost.seq.SequenceHandle`


  .. method:: Save(filename)

    Saves constraints to disk. So it will throw an error if you never called 
    the CalculateConstraints function on this DisCoContainer instance. 
    All added positions are lost in this process, only the actual constraints
    are saved.

    :param filename:    Name of file to be generated
    :type filename:     :class:`str`
    :raises:            :class:`RuntimeError` if *filename* does not exist
                        or when distance constraints have not been calculated 
                        yet


  .. method:: Load(filename)

    Loads a :class:`DisCoContainer` from disk and returns a readonly 
    :class:`DisCoContainer`. No positions can be added to this container 
    anymore, it just contains the actual constraints

    :param filename:    File to be loaded
    :type filename:     :class:`str` 
    :rtype:             :class:`DisCoContainer`
    :raises:            :class:`RuntimeError` if *filename* cannot be opened


  .. method:: AddData(aln, positions, pos_seqres_mapping)

    Adds template information. Call this function for every template you
    want to add.

    :param aln:         The alignment that will be used for sequence based 
                        clustering when calculating the constraints. The first
                        sequence must be the SEQRES given at initialization of 
                        the :class:`DisCoContainer`. The second sequence can be 
                        the ATOMSEQ or the SEQRES of the added template. This 
                        only influences the sequence based clustering.

    :param positions:   The CA positions of that template

    :param pos_seqres_mapping: A mapping to relate the *positions* to the internal
                               SEQRES. For every element in *positions* you must 
                               provide its residue number, i.e. the indexing in
                               the SEQRES starting from one.

    :type aln:          :class:`ost.seq.AlignmentHandle`
    :type positions:    :class:`ost.geom.Vec3List`
    :type pos_seqres_mapping: :class:`list` of :class:`int`

    :raises:            :class:`RuntimeError` if you already loaded this 
                        :class:`DisCoContainer` from disk (the thing is readonly)
                        or when an element in *pos_seqres* contains an invalid
                        residue number or when the first sequence in *aln* does
                        not exactly match the internal SEQRES (gaps are removed
                        for this check).


  .. method:: CalculateConstraints([dist_cutoff=15.0, gamma=70.0, \
                                    seqsim_clustering_cutoff=0.5, bin_size=0.5])

    Calculates the constraints using all the template information you already 
    added. The final constraints are lookup tables with equidistant bins as given
    by *bin_size*. They not only cover the range defined *dist_cutoff*, but add 
    3 Angstroms to allow the constraint to fade out.

    In principle constraint generation follows following steps (check out the paper
    for detailled information):

    * Generate one big multiple sequence alignment from all templates with the
      ost.seq.alg.MergePairwiseAlignments functionality using *seqres* as reference
    * Calculate all pairwise sequence similarities using BLOSUM62 (normalized =>
      (avg_seq_sim - min(BLOSUM62)) / (max(BLOSUM62) - min(BLOSUM62)))
    * Hierarchical clustering based on sequence similarities
    * Construct the distance constraints for every residue pair (*i*,*j*) as described
      in the paper, no constraint (*i*,*j*) will be generated if no pairwise distance 
      between (*i*,*j*) < *dist_cutoff* is observed in any template

    :param dist_cutoff: All pairwise distance constraints up to this cutoff will be
                        used.
    :param gamma:       This parameter influences the contribution of more remote 
                        homologues to the desired target sequence. The higher 
                        this value, the lower the contribution of more remote
                        homologues (where evolutionary distance is measured by
                        sequence similarity).
    :param seqsim_clustering_cutoff: All template sequences are hierarchically 
                        clustered by pairwise normalized sequence similarity 
                        (as given by BLOSUM62). The distance criterion is the 
                        average pairwise sequence similarity 
                        and clusters are merged only of this distance is below *seqsim_clustering_cutoff*
    :param bin_size:    The constraints are evaluated for equidistant bins of
                        of size *bin_size*.

    :type dist_cutoff:  :class:`float`
    :type gamma:        :class:`float`
    :type seqsim_clustering_cutoff:  :class:`float`
    :type bin_size:     :class:`float`


  .. method:: HasConstraint(i, j)

    Checks whether there is a constraint between the residues defined 
    by *i* and *j*. There can only be a constraint if you already called
    the CalculateConstraints function on this DisCoContainer instance and
    there has been at least one distance below  *dist_cutoff* in any of 
    the templates.

    :param i:           Residue number of the first residue (position in SEQRES
                        with indexing starting at 1)
    :param j:           Same for second residue

    :type i:            :class:`int`
    :type j:            :class:`int`
    :rtype:             :class:`bool`

    :raises:            :class:`RuntimeError` if *i* or *j* describe out of range 
                        positions

  .. method:: GetConstraint(i, j)

    Getter for the internal lookup table that describes the distance constraint
    between residues i and j.

    :param i:           Residue number of the first residue (position in SEQRES
                        with indexing starting at 1)
    :param j:           Same for second residue

    :type i:            :class:`int`
    :type j:            :class:`int`
    :rtype:             :class:`list` of :class:`float`

    :raises:            :class:`RuntimeError` if *i* or *j* describe out of range 
                        positions


  .. method:: GetScores(view)

    Evaluates the DisCo scores for all residues in the entity view. 
    The view must contain exactly one chain and all its residues must exactly 
    match the internal SEQRES with respect to residue number / one letter code.

    :param view:        View to evaluate
    :type view:         :class:`ost.mol.EntityView`

    :returns:           A list with one element per residue with the 
                        corresponding DisCo score
    :rtype:             :class:`list` of :class:`float`

    :raises:            :class:`RuntimeError` if there is more than one chain or
                        when any of the residues has a residue number / one 
                        letter that is inconsistent with the internal SEQRES.


  .. method:: GetScoresWithData(view)

    Evaluates the DisCo scores for all residues in the entity view. 
    The view must contain exactly one chain and all its residues must exactly 
    match the internal SEQRES with respect to residue number / one letter code.
    This is not everything... The final combination to QMEANDisCo requires more
    data from the DisCoContainer.


    :param view:        View to evaluate
    :type view:         :class:`ost.mol.EntityView`

    :returns:           A dictionary of lists, that contains not only the DisCo
                        scores but also other values. The lists can be accessed
                        by following keys:

                        * scores: The DisCo scores
                        * counts: The number of other residues within 15 A in 
                          the *view* with a valid constraint.
                        * avg_num_clusters: Average number of clusters observed
                          in all the constraints this residue is involved in.
                        * avg_max_seqsim: Average of highest per cluster 
                          sequence similarity observed in all the constraints 
                          this residue is involved in.
                        * avg_max_seqid: Average of highest per cluster 
                          sequence identity observed in all the constraints this 
                          residue is involved in.
                        * avg_variance: Average variance of the observed 
                          distances in all the constraints this residue is 
                          involved in.
                        * num_constraints: Total number of constraints that are 
                          present for this residue. Independent from the
                          inclusion radius that affects *counts*.


    :rtype:             :class:`dict` with :class:`str` as keys and 
                        :class:`list` of :class:`float` as values

    :raises:            :class:`RuntimeError` if there is more than one chain or
                        when any of the residues has a residue number / one 
                        letter that is inconsistent with the internal SEQRES.


  .. method:: GetSeqres()

    Getter for the SEQRES that has been set at initialization

    :rtype:             :class:`ost.seq.SequenceHandle`


  .. method:: GetDistCutoff()

    Getter for the *dist_cutoff* parameter you passed at the CalculateConstraints 
    function

    :rtype:             :class:`float`


  .. method:: GetGamma()

    Getter for the *gamma* parameter you passed at the CalculateConstraints 
    function

    :rtype:             :class:`float`


  .. method:: GetSeqSimClusteringCutoff()

    Getter for the *seqsim_clustering_cutoff* parameter you passed at the 
    CalculateConstraints function

    :rtype:             :class:`float`


  .. method:: GetBinSize()

    Getter for the *bin_size* parameter you passed at the CalculateConstraints 
    function

    :rtype:             :class:`float`
    

  .. method:: GetNumBins()

    The number of bins for the lookup table of each constraint. Its basically 
    the result of ceil((GetDistCutoff() + 3.0) / GetBinSize())

    :rtype:             :class:`int`

