Statistical Potentials of Mean Force
================================================================================

.. currentmodule:: qmean

Most of the available quality predictors are based on statistical potentials
of mean force as described by Sippl [sippl1990]_. 
The training and the application of the potentials of mean force 
is split into two objects. A statistic object is used to extract structural
features from training protein structures and internally stores them
using histogram techniques. A potential object then converts these
histograms into potential functions and assigns pseudo 
energies/scores to protein structures of interest. 
The number of available potentials has grown over the years. Only a subset 
of them is actually used in the official QMEAN scoring function.


* :ref:`interaction-potential`: 
                                Evaluates a pseudo energy
                                based on pairwise distances between all
                                chemically distinguishable heavy atoms.


* :ref:`cbeta-potential`: 
                          Evaluates a pseudo energy
                          based on pairwise distances between CB atoms.


* :ref:`reduced-potential`: 
                            Evaluates a pseudo energy between the
                            reduced representation of residues.
                            Every residue gets represented by its CA
                            position p and a directional component 
                            v = norm(ca_pos-n_pos) + norm(ca_pos-c_pos).
                            Assuming an interaction between residues r1 and r2, 
                            we can define a line l between p1 and p2. The 
                            potential then considers:

                            * dist => distance between p1 and p2
                            * alpha => angle between v1 and l
                            * beta => angle between v2 and l
                            * gamma => dihedral between (p1+v1,p1,p2,p2+v2)              

* :ref:`torsion-potential`: 
                            Evaluates a pseudo energy based on the identity of 
                            three consecutive residues and their phi/psi 
                            dihedral angles.


* :ref:`packing-potential`: 
                            Evaluates a pseudo energy by counting the number of 
                            surrounding atoms around all chemically 
                            distinguishable heavy atoms not belonging to the 
                            assessed residue itself.


* :ref:`cbpacking-potential`: 
                              Evaluates a pseudo energy by counting the number 
                              of other CB positions within a certain cutoff 
                              radius of the CB position of the residue to be 
                              evaluated.


* :ref:`hbond-potential`: 
                          Evaluates HBond pseudo energies similar to 
                          the one defined in the Rosetta energy 
                          function [kortemme2003]_. 
                          It considers the CA, C and O positions from backbone 
                          hbond acceptors in interaction with the N and H 
                          positions from the backbone hbond donors. 4 Parameters 
                          describe their relative orientation.

                          * dist => H-O distance
                          * alpha => O-H-N angle
                          * beta => C-N-H angle
                          * gamma => CA-C-O-H dihedral angle

                          A pseudo energy function for these parameters is 
                          evaluated for three different states. State 1 for 
                          helical residues, state 2 for extended residues and 
                          state 0 for other residues. If the state of two 
                          interacting particles is the same, thats the one from 
                          which the energy is extracted. In all other cases, 
                          the energy is extracted from the 0 state.

All required statistics and potentials are derived from the according base classes.
They are not accessible from Python.

.. class:: StatisticBase

.. class:: PotentialBase





.. _interaction-potential:

InteractionPotential
--------------------------------------------------------------------------------


.. literalinclude:: example_scripts/interaction_potential_example.py

.. class:: InteractionOpts

  Internal option class without init function. This class gets built when
  initializing the according statistic object and is automatically passed over
  when creating the potential. The class is considered to be readonly.

  .. attribute:: lower_cutoff

    Minimal distance two atom positions must have to be considered

  .. attribute:: upper_cutoff

    Maximal distance two atom positions can have to be considered

  .. attribute:: number_of_bins

    Number of equidistant bins building the internal histograms

  .. attribute:: sequence_sep

    Only atom Positions from residues separated by this amount of residues
    are considered

  .. attribute:: sigma

    Parameter to control the sparse data handling


.. class:: InteractionStatistic(l_cutoff, u_cutoff, nob, ssep)

  Builds histograms for each pair of all chemically distinguishable 
  atoms given the input parameters that can be filled using structural data.

  :param l_cutoff:      Minimal distance two atom positions must have to be 
                        considered 
  :param u_cutoff:      Maximal distance two atom positions can have to be 
                        considered
  :param nob:           Number of equidistant bins building the internal 
                        histograms
  :param ssep:          Only atom positions from residues separated by this 
                        amount of residues are considered

  :type l_cutoff:      :class:`float`
  :type u_cutoff:      :class:`float`
  :type nob:           :class:`int`
  :type ssep:          :class:`int`

  .. method:: Load(filename)

    Loads a :class:`InteractionStatistic` from disk

    :param filename:    Path to file containing the statistic
    :type filename:     :class:`str`

    :returns:           The loaded statistic
    :rtype:             :class:`InteractionStatistic`
    :raises:            :class:`RuntimeError` if **filename** does not exist

  .. method:: Save(filename)

    Saves a the statistic down to disk

    :param filename:    Path of the file, the statistic should be saved in
    :type filename:     :class:`str`

  .. method:: Extract(target, env, [weight = 1.0])

    Adds the value of the **weight** parameter to the according bins in the
    internal histograms

    :param target:      The target structure from which the statistics are 
                        extracted
    :param env:         The structure a particular atom position in the target
                        actually sees
    :param weight:      The value that gets added to the particular histogram 
                        positions

    :type target:       :class:`ost.mol.EntityView`
    :type env:          :class:`ost.mol.EntityView`
    :type weight:       :class:`float`

  .. method:: GetTotalCount()

    :returns:           The summed value in all distance bins in all 
                        pairwise histograms
    :rtype:             :class:`float`

  .. method:: GetCount(atom_one, atom_two, dist_bin)

    :param atom_one:    Identity of first atom
    :param atom_two:    Identity of second atom
    :param dist_bin:    Distance bin to be considered

    :type atom_one:     :ref:`ChemType <ChemType>`
    :type atom_two:     :ref:`ChemType <ChemType>`
    :type dist_bin:     :class:`int`

    :returns:           The histogram count for the particular pair of
                        atoms in the particular distance bin
    :rtype:             :class:`float`

  .. method:: GetCount(atom_one, atom_two)

    :param atom_one:    Identity of first atom
    :param atom_two:    Identity of second atom

    :type atom_one:     :ref:`ChemType <ChemType>`
    :type atom_two:     :ref:`ChemType <ChemType>`

    :returns:           The histogram count for the particular pair of
                        atoms
    :rtype:             :class:`float`

  .. method:: GetCount(dist_bin)

    :param dist_bin:    Distance bin to considered
    :type dist_bin:     :class:`int`

    :returns:           The histogram count for this particular
                        distance bins in all pairwise histograms

  .. method:: GetOptions()

    :returns:           A readonly option object
    :rtype:             :class:`InteractionOpts`


.. class:: InteractionPotential()

  Uses a :class:`InteractionStatistic` to generate a potential. The object must be
  created with its static Create function.

  .. method:: Load(filename)

    Loads a :class:`InteractionPotential` from disk

    :param filename:    Path to file containing the potential
    :type filename:     :class:`str`

    :returns:           The loaded potential
    :rtype:             :class:`InteractionPotential`
    :raises:            :class:`RuntimeError` if **filename** does not exist

  .. method:: Save(filename)

    Saves a the potential down to disk

    :param filename:    Path of the file, the potential should be saved in
    :type filename:     :class:`str`  

  .. method:: Create(stat, [sigma = 0.02, reference_state = "classic"])

    Uses the passed statistic to generate probability density functions, that
    are turned into a potential using the formalism described by Sippl.

    :param stat:        The statistic from which you want to create the potential
    :param sigma:       The sigma parameter to handle sparse data
    :param reference_state: The reference state to be used. Currently there
                            is only "classic" available, which implements
                            the reference state described by Sippl

    :type stat:         :class:`InteractionStatistic`
    :type sigma:        :class:`float`
    :type reference_state: :class:`str`

    :returns:           A newly created potential
    :rtype:             :class:`InteractionPotential`
    :raises:            :class:`RuntimeError` in case of invalid **reference_state**
    
      
  .. method:: Create(stat_one, stat_two, [sigma = 0.02, reference_state = "classic", max_energy=8.0)
    
    Uses the passed statistics to generate probability density functions, that
    are turned into a potential using the formalism developed for QMEANBrane.
    **stat_one** is considered to be saturated, whereas **stat_two** is 
    considered to be sparse.

    :param stat_one:    The saturates statistic from which you want to create
                        the potential
    :param stat_two:    The not necessarily saturated statistic from which you
                        want to create the potential
    
    :param sigma:       The sigma parameter to handle sparse data
    :param reference_state: The reference state to be used. Currently there
                            is only "classic" available, which implements
                            the reference state described by Sippl
    :param max_energy:  Due to the mathematical formalism, you will get infinit
                        scores for low pairwise distances
                        (Except your filling your stats with complete crap).
                        The max energy serves as a cap value.  

    :type stat_one:     :class:`InteractionStatistic`
    :type stat_two:     :class:`InteractionStatistic`
    :type sigma:        :class:`float`
    :type reference_state: :class:`str`
    :type max_energy:   :class:`float`

    :returns:           A newly created potential
    :rtype:             :class:`InteractionPotential`
    :raises:            :class:`RuntimeError` in case of invalid **reference_state**
                        or when the parametrization of the two stat objects is
                        inconsistent.

  .. method:: GetEnergy(atom_one, atom_two, dist)

    Get an energy value for a particular pair of atoms with their 
    positions at a certain distance

    :param atom_one:       First interaction partner
    :param atom_two:       Second interaction partner
    :param dist:           Distance of the two interaction partners positions

    :type amino_acid_one:  :ref:`ChemType <ChemType>`
    :type amino_acid_two:  :ref:`ChemType <ChemType>`
    :type dist:            :class:`float`

    :returns:              The energy, NaN if one of the partners is invalid or 
                           the distance is out of range
    :rtype:                :class:`float`

  .. method:: GetEnergy(res, env [normalize = True] )

    Gets the energy of **res** given the provided environment by evaluating
    all pairwise distances to the given environment.

    :param res:         Residue to be evaluated
    :param env:         The environment, the residue actually sees
    :param normalize:   Whether the final summed energy should be normalized by
                        total number of observed interactions
    :type res:          :class:`ost.mol.ResidueView`
    :type env:          :class:`ost.mol.EntityView`
    :type normalize:    :class:`bool`

    :returns:           The requested energy, NaN if no valid pairwise distance
                        is evaluated towards the environment
    :rtype:             :class:`float`

  .. method:: GetTotalEnergy(target, env, [normalize = True])

    Gets the energy of **target** given the provided environment by evaluating
    all pairwise distances from all residues of the **target** to the given 
    environment.

    :param target:      Structure to be evaluated
    :param env:         The environment, the residues from **target** actually sees
    :normalize:         Whether the final summed energy should be normalized by
                        total number of observed interactions
    :type target:       :class:`ost.mol.EntityView`
    :type env:          :class:`ost.mol.EntityView`
    :type normalize:    :class:`bool`

    :returns:           The requested energy, NaN if no valid pairwise distance
                        is evaluated towards the environment
    :rtype:             :class:`float`

  .. method:: GetEnergies(target, env, [normalize = True])

    Gets energies for every single residue in the **target** given the provided
    environment. For every residue, all pairwise interactions towards the
    environment are evaluated and summed up.

    :param target:      Structure to be evaluated
    :param env:         The environment, the residues from **target** actually sees
    :normalize:         Whether the final summed energy for every should be normalized by
                        total number of observed interactions for this residue
    :type target:       :class:`ost.mol.EntityView`
    :type env:          :class:`ost.mol.EntityView`
    :type normalize:    :class:`bool`

    :returns:           The requested energies, single elements can be NaN if no 
                        valid pairwise distance is evaluated towards the environment
                        for one particular residue
    :rtype:             :class:`list` of :class:`float`

  .. method:: GetOptions()

    :returns:           A readonly option object
    :rtype:             :class:`InteractionOpts`
     

.. _cbeta-potential:

CBetaPotential
--------------------------------------------------------------------------------

.. literalinclude:: example_scripts/cbeta_potential_example.py

.. class:: CBetaOpts

  Internal option class without init function. This class gets built when
  initializing the according statistic object and is automatically passed over
  when creating the potential. The class is considered to be readonly.

  .. attribute:: lower_cutoff

    Minimal distance two CB positions must have to be considered

  .. attribute:: upper_cutoff

    Maximal distance two CB positions can have to be considered

  .. attribute:: number_of_bins

    Number of equidistant bins building the internal histograms

  .. attribute:: sequence_sep

    Only CB Positions from residues separated by this amount of residues
    are considered

  .. attribute:: sigma

    Parameter to control the sparse data handling


.. class:: CBetaStatistic(l_cutoff, u_cutoff, nob, ssep)

  Builds histograms for each pair of the 20 standard amino acids given the 
  input parameters that can be filled using structural data.

  :param l_cutoff:      Minimal distance two CB positions must have to be 
                        considered 
  :param u_cutoff:      Maximal distance two CB positions can have to be 
                        considered
  :param nob:           Number of equidistant bins building the internal 
                        histograms
  :param ssep:          Only CB Positions from residues separated by this amount of residues
                        are considered

  :type l_cutoff:      :class:`float`
  :type u_cutoff:      :class:`float`
  :type nob:           :class:`int`
  :type ssep:          :class:`int`

  .. method:: Load(filename)

    Loads a :class:`CBetaStatistic` from disk

    :param filename:    Path to file containing the statistic
    :type filename:     :class:`str`

    :returns:           The loaded statistic
    :rtype:             :class:`CBetaStatistic`
    :raises:            :class:`RuntimeError` if **filename** does not exist

  .. method:: Save(filename)

    Saves a the statistic down to disk

    :param filename:    Path of the file, the statistic should be saved in
    :type filename:     :class:`str`

  .. method:: Extract(target, env, [weight = 1.0])

    Adds the value of the **weight** parameter to the according bins in the
    internal histograms

    :param target:      The target structure from which the statistics are 
                        extracted
    :param env:         The structure a particular CB Position in the target
                        actually sees
    :param weight:      The value that gets added to the particular histogram 
                        positions

    :type target:       :class:`ost.mol.EntityView`
    :type env:          :class:`ost.mol.EntityView`
    :type weight:       :class:`float`

  .. method:: GetTotalCount()

    Sums up all values in all distance bins in all pairwise histograms

    :returns:           The summed histogram value
    :rtype:             :class:`int`

  .. method:: GetCount(aa_one, aa_two, dist_bin)

    :param aa_one:      Identity of first amino acid
    :param aa_two:      Identity of second amino acid
    :param dist_bin:    Distance bin to be considered

    :type aa_one:       :class:`ost.conop.AminoAcid`
    :type aa_two:       :class:`ost.conop.AminoAcid`
    :type dist_bin:     :class:`int`

    :returns:           The histogram count for the particular pair of
                        amino acids in the particular distance bin
    :rtype:             :class:`int`

  .. method:: GetCount(aa_one, aa_two)

    :param aa_one:      Identity of first amino acid
    :param aa_two:      Identity of second amino acid

    :type aa_one:       :class:`ost.conop.AminoAcid`
    :type aa_two:       :class:`ost.conop.AminoAcid`

    :returns:           The histogram count for the particular pair of
                        amino acids
    :rtype:             :class:`int`


  .. method:: GetCount(dist_bin)

    :param dist_bin:    Distance bin to considered
    :type dist_bin:     :class:`int`

    :returns:          The histogram count for this particular distance bin
                       summed over all pairwise histograms

  .. method:: GetOptions()

    :returns:           A readonly option object
    :rtype:             :class:`CBetaOpts`


.. class:: CBetaPotential()

  Uses a :class:`CBetaStatistic` to generate a potential. The object must be
  created with its static Create function.

  .. method:: Load(filename)

    Loads a :class:`CBetaPotential` from disk

    :param filename:    Path to file containing the potential
    :type filename:     :class:`str`

    :returns:           The loaded potential
    :rtype:             :class:`CBetaPotential`
    :raises:            :class:`RuntimeError` if **filename** does not exist

  .. method:: Save(filename)

    Saves a the potential down to disk

    :param filename:    Path of the file, the potential should be saved in
    :type filename:     :class:`str`  

  .. method:: Create(stat, [sigma = 0.02, reference_state = "classic"])

    Uses the passed statistic to generate probability density functions, that
    are turned into a potential using the formalism described by Sippl.

    :param stat:        The statistic from which you want to create the potential
    :param sigma:       The sigma parameter to handle sparse data
    :param reference_state: The reference state to be used. Currently there
                            is only "classic" available, which implements
                            the reference state described by Sippl

    :type stat:         :class:`CBetaStatistic`
    :type sigma:        :class:`float`
    :type reference_state: :class:`str`

    :returns:           A newly created potential
    :rtype:             :class:`CBetaPotential`
    :raises:            :class:`RuntimeError` in case of invalid **reference_state**
    
      
  .. method:: Create(stat_one, stat_two, [sigma = 0.02, reference_state = "classic", max_energy=8.0)
    
    Uses the passed statistics to generate probability density functions, that
    are turned into a potential using the formalism developed for QMEANBrane.
    **stat_one** is considered to be saturated, whereas **stat_two** is 
    considered to be sparse.

    :param stat_one:    The saturates statistic from which you want to create
                        the potential
    :param stat_two:    The not necessarily saturated statistic from which you
                        want to create the potential
    
    :param sigma:       The sigma parameter to handle sparse data
    :param reference_state: The reference state to be used. Currently there
                            is only "classic" available, which implements
                            the reference state described by Sippl
    :param max_energy:  Due to the mathematical formalism, you will get infinit
                        scores for low pairwise distances
                        (Except your filling your stats with complete crap).
                        The max energy serves as a cap value.  

    :type stat_one:     :class:`CBetaStatistic`
    :type stat_two:     :class:`CBetaStatistic`
    :type sigma:        :class:`float`
    :type reference_state: :class:`str`
    :type max_energy:   :class:`float`

    :returns:           A newly created potential
    :rtype:             :class:`CBetaPotential`
    :raises:            :class:`RuntimeError` in case of invalid **reference_state**
                        or when the parametrization of the two stat objects is
                        inconsistent.

  .. method:: GetEnergy(amino_acid_one, amino_acid_two, dist)

    Get an energy value for a particular pair of amino acids with their 
    CB positions at a certain distance

    :param amino_acid_one: First interaction partner
    :param amino_acid_two: Second interaction partner
    :param dist:           Distance of the two interaction partners CB positions

    :type amino_acid_one:  :class:`ost.conop.AminoAcid`
    :type amino_acid_two:  :class:`ost.conop.AminoAcid`
    :type dist:            :class:`float`

    :returns:              The energy, NaN if one of the partners is invalid or 
                           the distance is out of range
    :rtype:                :class:`float`

  .. method:: GetEnergy(res, env [normalize = True] )

    Gets the energy of **res** given the provided environment by evaluating
    all pairwise distances to the given environment.

    :param res:         Residue to be evaluated
    :param env:         The environment, the residue actually sees
    :param normalize:   Whether the final summed energy should be normalized by
                        total number of observed interactions
    :type res:          :class:`ost.mol.ResidueView`
    :type env:          :class:`ost.mol.EntityView`
    :type normalize:    :class:`bool`

    :returns:           The requested energy, NaN if no valid pairwise distance
                        is evaluated towards the environment
    :rtype:             :class:`float`

  .. method:: GetTotalEnergy(target, env, [normalize = True])

    Gets the energy of **target** given the provided environment by evaluating
    all pairwise distances from all residues of the **target** to the given 
    environment.

    :param target:      Structure to be evaluated
    :param env:         The environment, the residues from **target** actually sees
    :normalize:         Whether the final summed energy should be normalized by
                        total number of observed interactions
    :type target:       :class:`ost.mol.EntityView`
    :type env:          :class:`ost.mol.EntityView`
    :type normalize:    :class:`bool`

    :returns:           The requested energy, NaN if no valid pairwise distance
                        is evaluated towards the environment
    :rtype:             :class:`float`

  .. method:: GetEnergies(target, env, [normalize = True])

    Gets energies for every single residue in the **target** given the provided
    environment. For every residue, all pairwise interactions towards the
    environment are evaluated and summed up.

    :param target:      Structure to be evaluated
    :param env:         The environment, the residues from **target** actually sees
    :normalize:         Whether the final summed energy for every should be normalized by
                        total number of observed interactions for this residue
    :type target:       :class:`ost.mol.EntityView`
    :type env:          :class:`ost.mol.EntityView`
    :type normalize:    :class:`bool`

    :returns:           The requested energies, single elements can be NaN if no 
                        valid pairwise distance is evaluated towards the environment
                        for one particular residue
    :rtype:             :class:`list` of :class:`float`

  .. method:: GetOptions()

    :returns:           A readonly option object
    :rtype:             :class:`CBetaOpts`



.. _reduced-potential:

ReducedPotential
--------------------------------------------------------------------------------


.. class:: ReducedOpts

  Internal option class without init function. This class gets built when
  initializing the according statistic object and is automatically passed over
  when creating the potential. The class is considered to be readonly.

  .. attribute:: lower_cutoff

    Minimal distance two CA positions must have to be considered

  .. attribute:: upper_cutoff

    Maximal distance two CA positions can have to be considered

  .. attribute:: num_angle_bins

    Number of equidistant angle bins building the internal histograms

  .. attribute:: num_dihedral_bins

    Number of equidistance dihedral bins building the internal histograms

  .. attribute: num_dist_bins

    Number of equidistance distance bins building the internal histograms

  .. attribute:: sequence_sep

    Only pairwise interactions from residues separated by this amount of 
    residues are considered

  .. attribute:: sigma

    Parameter to control the sparse data handling


.. class:: ReducedStatistic(l_cutoff, u_cutoff, nab, ndab, ndb, ssep)

  Builds histograms for each pair of the 20 standard amino acids given the 
  input parameters that can be filled using structural data.

  :param l_cutoff:      Minimal distance two CA positions must have to be 
                        considered 
  :param u_cutoff:      Maximal distance two CA positions can have to be 
                        considered
  :param nab:           Number of equidistant angle bins building the internal 
                        histograms
  :param ndab:          Number of equidistant dihedral angle bins building the internal 
                        histograms
  :param ndb:           Number of equidistant distance bins building the internal 
                        histograms
  :param ssep:          Only CA positions from residues separated by this amount of residues
                        are considered

  :type l_cutoff:      :class:`float`
  :type u_cutoff:      :class:`float`
  :type nab:           :class:`int`
  :type ndab:          :class:`int`
  :type ndb:           :class:`int`
  :type ssep:          :class:`int`

  .. method:: Load(filename)

    Loads a :class:`ReducedStatistic` from disk

    :param filename:    Path to file containing the statistic
    :type filename:     :class:`str`

    :returns:           The loaded statistic
    :rtype:             :class:`ReducedStatistic`
    :raises:            :class:`RuntimeError` if **filename** does not exist

  .. method:: Save(filename)

    Saves the statistic down to disk

    :param filename:    Path of the file, the statistic should be saved in
    :type filename:     :class:`str`

  .. method:: Extract(target, env, [weight = 1.0])

    Adds the value of the **weight** parameter to the according bins in the
    internal histograms

    :param target:      The target structure from which the statistics are 
                        extracted
    :param env:         The structure a particular residue in the target
                        actually sees
    :param weight:      The value that gets added to the particular histogram 
                        positions

    :type target:       :class:`ost.mol.EntityView`
    :type env:          :class:`ost.mol.EntityView`
    :type weight:       :class:`float`

  .. method:: GetTotalCount()

    Sums up all values in all bins of the histograms

    :returns:           The summed histogram value
    :rtype:             :class:`float`

  .. method:: GetCount(aa_one, aa_two, dist_bin, alpha_bin, beta_bin, gamma_bin)

    :param aa_one:      Identity of first amino acid
    :param aa_two:      Identity of second amino acid
    :param dist_bin:    Distance bin to be considered
    :param alpha_bin:   Alpha angle bin to be considered
    :param beta_bin:    Beta angle bin to be considered
    :param gamma_bin:   Gamma Dihedral angle bin to be considered

    :type aa_one:       :class:`ost.conop.AminoAcid`
    :type aa_two:       :class:`ost.conop.AminoAcid`
    :type dist_bin:     :class:`int`
    :type alpha_bin:    :class:`int`
    :type beta_bin:     :class:`int`
    :type gamma_bin:    :class:`int`

    :returns:           The histogram count for the particular pair of
                        amino acids in the particular distance, angle and 
                        dihedral bin
    :rtype:             :class:`float`

  .. method:: GetCount(aa_one, aa_two)

    :param aa_one:      Identity of first amino acid
    :param aa_two:      Identity of second amino acid

    :type aa_one:       :class:`ost.conop.AminoAcid`
    :type aa_two:       :class:`ost.conop.AminoAcid`

    :returns:           The histogram count for the particular pair of
                        amino acids
    :rtype:             :class:`float`


  .. method:: GetCount(dist_bin, alpha_bin, beta_bin, gamma_bin)

    Sums up all values for one particular distance, angle and dihedral bin 
    for all possible pairs of amino acids

    :param dist_bin:    Distance bin to considered
    :type dist_bin:     :class:`int`

  .. method:: GetOptions()

    :returns:           A readonly option object
    :rtype:             :class:`ReducedOpts`


.. class:: ReducedPotential()

  Uses a :class:`ReducedStatistic` to generate a potential. The object must be
  created with its static Create function.

  .. method:: Load(filename)

    Loads a :class:`ReducedPotential` from disk

    :param filename:    Path to file containing the potential
    :type filename:     :class:`str`

    :returns:           The loaded potential
    :rtype:             :class:`ReducedPotential`
    :raises:            :class:`RuntimeError` if **filename** does not exist

  .. method:: Save(filename)

    Saves a the potential down to disk

    :param filename:    Path of the file, the potential should be saved in
    :type filename:     :class:`str`  

  .. method:: Create(stat, [sigma = 0.02, reference_state = "classic"])

    Uses the passed statistic to generate probability density functions, that
    are turned into a potential using the formalism described by Sippl.

    :param stat:        The statistic from which you want to create the potential
    :param sigma:       The sigma parameter to handle sparse data
    :param reference_state: The reference state to be used. Currently there
                            is only "classic" available, which implements
                            the reference state described by Sippl

    :type stat:         :class:`ReducedStatistic`
    :type sigma:        :class:`float`
    :type reference_state: :class:`str`

    :returns:           A newly created potential
    :rtype:             :class:`ReducedPotential`
    :raises:            :class:`RuntimeError` in case of invalid **reference_state**
    

  .. method:: GetEnergy(amino_acid_one, amino_acid_two, dist, alpha, beta, gamma)

    Get an energy value for a particular pair of amino acids with their 
    their relative orientation described by dist, alpha beta and gamma

    :param amino_acid_one: First interaction partner
    :param amino_acid_two: Second interaction partner
    :param dist:           Distance of the two interaction partners CA positions
    :param alpha:          First angle describing the relative orientation
    :param beta:           First angle describing the relative orientation
    :param gamma:          Dihedral angle describing the relative orientation

    :type amino_acid_one:  :class:`ost.conop.AminoAcid`
    :type amino_acid_two:  :class:`ost.conop.AminoAcid`
    :type dist:            :class:`float`
    :type alpha:           :class:`float`
    :type beta:            :class:`float`
    :type gamma:           :class:`float`

    :returns:              The energy, NaN if one of the partners is invalid or 
                           the distance is out of range
    :rtype:                :class:`float`

  .. method:: GetEnergy(res, env [normalize = True] )

    Gets the energy of **res** given the provided environment by evaluating
    all pairwise distances to the given environment.

    :param res:         Residue to be evaluated
    :param env:         The environment, the residue actually sees
    :param normalize:   Whether the final summed energy should be normalized by
                        total number of observed interactions
    :type res:          :class:`ost.mol.ResidueView`
    :type env:          :class:`ost.mol.EntityView`
    :type normalize:    :class:`bool`

    :returns:           The requested energy, NaN if no valid pairwise distance
                        is evaluated towards the environment
    :rtype:             :class:`float`

  .. method:: GetTotalEnergy(target, env, [normalize = True])

    Gets the energy of **target** given the provided environment by evaluating
    all pairwise distances from all residues of the **target** to the given 
    environment.

    :param target:      Structure to be evaluated
    :param env:         The environment, the residues from **target** actually sees
    :normalize:         Whether the final summed energy should be normalized by
                        total number of observed interactions
    :type target:       :class:`ost.mol.EntityView`
    :type env:          :class:`ost.mol.EntityView`
    :type normalize:    :class:`bool`

    :returns:           The requested energy, NaN if no valid pairwise distance
                        is evaluated towards the environment
    :rtype:             :class:`float`

  .. method:: GetEnergies(target, env, [normalize = True])

    Gets energies for every single residue in the **target** given the provided
    environment. For every residue, all pairwise interactions towards the
    environment are evaluated and summed up.

    :param target:      Structure to be evaluated
    :param env:         The environment, the residues from **target** actually sees
    :normalize:         Whether the final summed energy for every should be normalized by
                        total number of observed interactions for this residue
    :type target:       :class:`ost.mol.EntityView`
    :type env:          :class:`ost.mol.EntityView`
    :type normalize:    :class:`bool`

    :returns:           The requested energies, single elements can be NaN if no 
                        valid pairwise distance is evaluated towards the environment
                        for one particular residue
    :rtype:             :class:`list` of :class:`float`

  .. method:: GetOptions()

    :returns:           A readonly option object
    :rtype:             :class:`ReducedOpts`


.. _torsion-potential:
   
TorsionPotential
--------------------------------------------------------------------------------

Since the torsion sampler considers triplets of amino acids, we need to define 
them. This is done with the so called torsion group definitions. Three strings 
represent the according positions of the consecutive amino acids. They are 
combined by “-”. It is either possible to use the keyword “all”, or write out 
all allowed amino acids by their three letter code and separate them by ”,”. 
An example would be: “all-VAL,ILE-PRO”. There are cases where a tripeptide can 
match several group definitions. The group definitions are stored internally as 
a list. This list is iterated at every evaluation of three consecutive amino 
acids and the first hit is decisive.



.. literalinclude:: example_scripts/torsion_potential_example.py


.. class:: TorsionOpts

  Internal option class without init function. This class gets built when
  initializing the according statistic object and is automatically passed over
  when creating the potential. The class is considered to be readonly.

  .. attribute:: number_of_bins

    List of length 6 describing the number of bins for each phi/psi dihedral
    angle for three consecutive residues

  .. attribute:: group_identifier

    List of definitions to group triplets of amino acids together.

  .. attribute:: sigma

    Parameter to control the sparse data handling


.. class:: TorsionStatistic(group_identifier, number_of_bins)

  Builds histograms for each group of amino acids given the 
  the **number_of_bins** that can be filled using structural data.

  :param group_identifier: Definition of amino acids that should be grouped 
                           together

  :param number_of_bins: Number of bins for all phi/psi dihedrals when looking
                         at three consecutive amino acids  

  :type group_identifier: :class:`list` of :class:`str`
  :type number_of_bins:   :class:`list` of 6 :class:`int` values

  .. method:: Load(filename)

    Loads a :class:`TorsionStatistic` from disk

    :param filename:    Path to file containing the statistic
    :type filename:     :class:`str`

    :returns:           The loaded statistic
    :rtype:             :class:`TorsionStatistic`
    :raises:            :class:`RuntimeError` if **filename** does not exist

  .. method:: Save(filename)

    Saves a the statistic down to disk

    :param filename:    Path of the file, the statistic should be saved in
    :type filename:     :class:`str`

  .. method:: Extract(target, [weight = 1.0])

    Adds the value of the **weight** parameter to the according bins in the
    internal histograms

    :param target:      The target structure from which the statistics are 
                        extracted
    :param weight:      The value that gets added to the particular histogram 
                        positions

    :type target:       :class:`ost.mol.EntityView`
    :type weight:       :class:`float`

  .. method:: GetTotalCount()

    Sums up all values in all histograms

    :returns:           The summed histogram value
    :rtype:             :class:`float`

  .. method:: GetTotalCount(group_identifier)

    Sums up all values in the histogram specific for **group_identifier**

    :param group_identifier: The group defining the amino acid triplet you're
                             interested in

    :type group_identifier: :class:`str`

    :returns:           The summed histogram value
    :rtype:             :class:`float`

  .. method:: GetCount(group_identifier, bins)

    :param group_identifier: The group defining the amino acid triplet you're
                             interested in       
    :param bins:        The phi/psi dihedral angle bins you're interested in

    :type group_identifier: :class:`list` of :class:`str`
    :type bins:         :class:`list` of 6 :class:`int` values

    :returns:           The histogram count for the particular group of
                        amino acid  triplets in the according dihedral
                        angle bins
    :rtype:             :class:`float`

  .. method:: GetCount(bins)

    :param bins:        The phi/psi bins you're interested in

    :type bins:         :class:`list` of 6 :class:`int` values         

    :returns:           The histogram count for the particular bin combination
                        summed over all possible groups of amino acid triplets
    :rtype:             :class:`float`


.. class:: TorsionPotential()

  Uses a :class:`TorsionStatistic` to generate a potential. The object must be
  created with its static Create function.

  .. method:: Load(filename)

    Loads a :class:`TorsionPotential` from disk

    :param filename:    Path to file containing the potential
    :type filename:     :class:`str`

    :returns:           The loaded potential
    :rtype:             :class:`TorsionPotential`
    :raises:            :class:`RuntimeError` if **filename** does not exist

  .. method:: Save(filename)

    Saves a the potential down to disk

    :param filename:    Path of the file, the potential should be saved in
    :type filename:     :class:`str`  

  .. method:: Create(stat, [sigma = 0.02, reference_state = "classic"])

    Uses the passed statistic to generate probability density functions, that
    are turned into a potential using the formalism described by Sippl. There
    will be created a separate potential for all groups.

    :param stat:        The statistic from which you want to create the potential
    :param sigma:       The sigma parameter to handle sparse data
    :param reference_state: The reference state to be used. Currently there
                            is "classic" available, which implements
                            the reference state described by Sippl, an
                            alternative is "uniform".

    :type stat:         :class:`TorsionStatistic`
    :type sigma:        :class:`float`
    :type reference_state: :class:`str`

    :returns:           A newly created potential
    :rtype:             :class:`TorsionPotential`
    :raises:            :class:`RuntimeError` in case of invalid **reference_state**
    
      
  .. method:: GetEnergy(group_identifier, angles)

    Get an energy value for a particular group of amino acid triplets given certain
    phi/psi dihedral angles 

    :param group_identifier: The group of amino acid triplets you're interested in
    :param angles:      The angles for the phi/psi backbone dihedral angles

    :type group_identifier: :class:`str`
    :type angles:       :class:`list` of 6 :class:`float` values

    :returns:              The energy
    :rtype:                :class:`float`

  .. method:: GetEnergy(res)

    Gets the energy of **res** 

    :param res:         Residue to be evaluated
    :type res:          :class:`ost.mol.ResidueView`

    :returns:           The requested energy, NaN if one or both torsions in the
                        residue are invalid. This happend if the residue or one
                        of the two neighbours is not connected to another 
                        residue
    :rtype:             :class:`float`

  .. method:: GetTotalEnergy(target, [normalize = True])

    Gets the energy of **target** by summing up all single residue scores

    :param target:      Structure to be evaluated
    :normalize:         Whether the final summed energy should be normalized by
                        total number of observed interactions
    :type target:       :class:`ost.mol.EntityView`
    :type normalize:    :class:`bool`

    :returns:           The requested energy, NaN if no valid interaction
                        is evaluated
    :rtype:             :class:`float`

  .. method:: GetEnergies(target)

    Gets energies for every single residue in the **target**. There might be NaN values
    for if there are invalid torsions involved in some residues (termini)

    :param target:      Structure to be evaluated
    :type target:       :class:`ost.mol.EntityView`

    :returns:           The requested energies
    :rtype:             :class:`list` of :class:`float`

  .. method:: GetOptions()

    :returns:           A readonly option object
    :rtype:             :class:`TorsionOpts`



.. _packing-potential:

PackingPotential
--------------------------------------------------------------------------------

.. literalinclude:: example_scripts/packing_potential_example.py

.. class:: PackingOpts

  Internal option class without init function. This class gets built when
  initializing the according statistic object and is automatically passed over
  when creating the potential. The class is considered to be readonly.

  .. attribute:: cutoff

    All atom positions within this radius around the atom position of interest
    are counted

  .. attribute:: max_counts

    The max amount of counts => cap value

  .. attribute:: sigma

    Parameter to control the sparse data handling


.. class:: PackingStatistic(cutoff, max_counts)

  Builds histograms for each pair of the 20 standard amino acids given the 
  input parameters that can be filled using structural data.

  :param cutoff:        Cutoff radius to count other atom positions in the 
                        environment
  :param max_counts:    Upper cap of other atom positions to be counted

  :type cutoff:         :class:`float`
  :type max_counts:     :class:`int`

  .. method:: Load(filename)

    Loads a :class:`PackingStatistic` from disk

    :param filename:    Path to file containing the statistic
    :type filename:     :class:`str`

    :returns:           The loaded statistic
    :rtype:             :class:`PackingStatistic`
    :raises:            :class:`RuntimeError` if **filename** does not exist

  .. method:: Save(filename)

    Saves the statistic down to disk

    :param filename:    Path of the file, the statistic should be saved in
    :type filename:     :class:`str`

  .. method:: Extract(target, env, [weight = 1.0])

    Adds the value of the **weight** parameter to the according bins in the
    internal histograms

    :param target:      The target structure from which the statistics are 
                        extracted
    :param env:         The structure a particular atom position in the target
                        actually sees
    :param weight:      The value that gets added to the particular histogram 
                        positions

    :type target:       :class:`ost.mol.EntityView`
    :type env:          :class:`ost.mol.EntityView`
    :type weight:       :class:`float`

  .. method:: GetTotalCount()

    Sums up all values in all histogram bins

    :returns:           The summed histogram value
    :rtype:             :class:`float`

  .. method:: GetCount(atom)

    :param atom:        Identity of the atom

    :type aa:           :ref:`ChemType <ChemType>`

    :returns:           The sum of all counts in the histogram for one 
                        particular atom

    :rtype:             :class:`float`

  .. method:: GetCount(count)

    :param count:       The particular count you're interested in

    :type count:       :class:`int`

    :returns:           The sum of all counts in the histogram for
                        this particular count

    :rtype:             :class:`float`


  .. method:: GetCount(atom, count)

    :param atom:        Identity of the atom
    :param count:       The particular count you're interested in

    :type atom:         :ref:`ChemType <ChemType>`
    :type count:        :class:`int`

    :returns:           The histogram value for the particular 
                        atom / count combination

    :rtype:             :class:`float`


  .. method:: GetOptions()

    :returns:           A readonly option object
    :rtype:             :class:`PackingOpts`



.. class:: PackingPotential()

  Uses a :class:`PackingStatistic` to generate a potential. The object must be
  created with its static Create function.

  .. method:: Load(filename)

    Loads a :class:`PackingPotential` from disk

    :param filename:    Path to file containing the potential
    :type filename:     :class:`str`

    :returns:           The loaded potential
    :rtype:             :class:`PackingPotential`
    :raises:            :class:`RuntimeError` if **filename** does not exist

  .. method:: Save(filename)

    Saves a the potential down to disk

    :param filename:    Path of the file, the potential should be saved in
    :type filename:     :class:`str`  

  .. method:: Create(stat, [sigma = 0.02, reference_state = "classic"])

    Uses the passed statistic to generate probability density functions, that
    are turned into a potential using the formalism described by Sippl.

    :param stat:        The statistic from which you want to create the potential
    :param sigma:       The sigma parameter to handle sparse data
    :param reference_state: The reference state to be used. Currently there
                            is "classic" available, which implements
                            the reference state described by Sippl, an
                            alternative is "uniform".

    :type stat:         :class:`PackingStatistic`
    :type sigma:        :class:`float`
    :type reference_state: :class:`str`

    :returns:           A newly created potential
    :rtype:             :class:`PackingPotential`
    :raises:            :class:`RuntimeError` in case of invalid **reference_state**
    

  .. method:: GetEnergy(atom, count)

    Get an energy value for a particular amino acids with a specific count

    :param atom:           The atom of interest
    :param count:          The count

    :type aa:              :ref:`ChemType <ChemType>`
    :type count:           :class:`float`

    :returns:              The energy given the input parameters
    :rtype:                :class:`float`

  .. method:: GetEnergy(res, env, normalize)

    Gets the energy of **res** given the provided environment.

    :param res:         Residue to be evaluated
    :param env:         The environment, the residue actually sees
    :normalize:         Whether the energy value should be normalized by the
                        number of atoms in the residue contributing to the final
                        energy
    :type res:          :class:`ost.mol.ResidueView`
    :type env:          :class:`ost.mol.EntityView`
    :type normalize:    :class:`bool`

    :returns:           The requested energy
    :rtype:             :class:`float`

  .. method:: GetTotalEnergy(target, env, normalize)

    Gets the energy of **target** given the provided environment

    :param target:      Structure to be evaluated
    :param env:         The environment, the residues from **target** actually sees
    :normalize:         Whether the summed energy should be normalized by the
                        total number of atoms contributing to the final
                        energy

    :type target:       :class:`ost.mol.EntityView`
    :type env:          :class:`ost.mol.EntityView`
    :type normalize:    :class:`bool`

    :returns:           The requested energy, NaN if no interaction could be observed
    :rtype:             :class:`float`

  .. method:: GetEnergies(target, env, normalize)

    Gets energies for every single residue in the **target** given the provided
    environment. 

    :param target:      Structure to be evaluated
    :param env:         The environment, the residues from **target** actually sees
    :param normalize:   Whether the per residue energy values should be 
                        normalized by the number of atoms in a residue contributing
                        to its final energy

    :type target:       :class:`ost.mol.EntityView`
    :type env:          :class:`ost.mol.EntityView`
    :type normalize:    :class:`bool`

    :returns:           The requested energies
    :rtype:             :class:`list` of :class:`float`

  .. method:: GetOptions()

    :returns:           A readonly option object
    :rtype:             :class:`PackingOpts`



.. _cbpacking-potential:

CBPackingPotential
--------------------------------------------------------------------------------


.. class:: CBPackingOpts

  Internal option class without init function. This class gets built when
  initializing the according statistic object and is automatically passed over
  when creating the potential. The class is considered to be readonly.

  .. attribute:: cutoff

    All CB positions within this radius around the CB position of interest
    are counted

  .. attribute:: max_counts

    The max amount of counts => cap value

  .. attribute:: sigma

    Parameter to control the sparse data handling


.. class:: CBPackingStatistic(cutoff, max_counts)

  Builds histograms for each of the 20 standard amino acids given the 
  input parameters that can be filled using structural data.

  :param cutoff:        Cutoff radius to count other CB positions in the 
                        environment
  :param max_counts:    Upper cap of other CB positions to be counted

  :type cutoff:         :class:`float`
  :type max_counts:     :class:`int`

  .. method:: Load(filename)

    Loads a :class:`CBPackingStatistic` from disk

    :param filename:    Path to file containing the statistic
    :type filename:     :class:`str`

    :returns:           The loaded statistic
    :rtype:             :class:`CBPackingStatistic`
    :raises:            :class:`RuntimeError` if **filename** does not exist

  .. method:: Save(filename)

    Saves the statistic down to disk

    :param filename:    Path of the file, the statistic should be saved in
    :type filename:     :class:`str`

  .. method:: Extract(target, env, [weight = 1.0])

    Adds the value of the **weight** parameter to the according bins in the
    internal histograms

    :param target:      The target structure from which the statistics are 
                        extracted
    :param env:         The structure a particular CB position in the target
                        actually sees
    :param weight:      The value that gets added to the particular histogram 
                        positions

    :type target:       :class:`ost.mol.EntityView`
    :type env:          :class:`ost.mol.EntityView`
    :type weight:       :class:`float`

  .. method:: GetTotalCount()

    Sums up all values in all histogram bins

    :returns:           The summed histogram value
    :rtype:             :class:`float`

  .. method:: GetCount(aa)

    :param aa:          Identity of the amino acid

    :type aa:           :class:`ost.conop.AminoAcid`

    :returns:           The sum of all counts in the histogram for one 
                        particular amino acid

    :rtype:             :class:`float`

  .. method:: GetCount(count)

    :param count:       The particular count you're interested in

    :type count:       :class:`int`

    :returns:           The sum of all counts in the histogram for
                        this particular count

    :rtype:             :class:`float`


  .. method:: GetCount(aa, count)

    :param aa:          Identity of the amino acid
    :param count:       The particular count you're interested in

    :type aa:           :class:`ost.conop.AminoAcid`
    :type count:        :class:`int`

    :returns:           The histogram value for the particular 
                        amino acid / count combination

    :rtype:             :class:`float`


  .. method:: GetOptions()

    :returns:           A readonly option object
    :rtype:             :class:`CBetaOpts`





.. class:: CBPackingPotential()

  Uses a :class:`CBPackingStatistic` to generate a potential. The object must be
  created with its static Create function.

  .. method:: Load(filename)

    Loads a :class:`CBPackingPotential` from disk

    :param filename:    Path to file containing the potential
    :type filename:     :class:`str`

    :returns:           The loaded potential
    :rtype:             :class:`CBPackingPotential`
    :raises:            :class:`RuntimeError` if **filename** does not exist

  .. method:: Save(filename)

    Saves a the potential down to disk

    :param filename:    Path of the file, the potential should be saved in
    :type filename:     :class:`str`  

  .. method:: Create(stat, [sigma = 0.02, reference_state = "classic"])

    Uses the passed statistic to generate probability density functions, that
    are turned into a potential using the formalism described by Sippl.

    :param stat:        The statistic from which you want to create the potential
    :param sigma:       The sigma parameter to handle sparse data
    :param reference_state: The reference state to be used. Currently there
                            is "classic" available, which implements
                            the reference state described by Sippl, an
                            alternative is "uniform".

    :type stat:         :class:`CBPackingStatistic`
    :type sigma:        :class:`float`
    :type reference_state: :class:`str`

    :returns:           A newly created potential
    :rtype:             :class:`CBPackingPotential`
    :raises:            :class:`RuntimeError` in case of invalid **reference_state**
    

  .. method:: GetEnergy(aa, count)

    Get an energy value for a particular amino acids with a specific count

    :param aa:             The amino acid of interest
    :param count:          The count

    :type aa:              :class:`ost.conop.AminoAcid`
    :type count:           :class:`float`

    :returns:              The energy given the input parameters
    :rtype:                :class:`float`

  .. method:: GetEnergy(res, env)

    Gets the energy of **res** given the provided environment.

    :param res:         Residue to be evaluated
    :param env:         The environment, the residue actually sees
    :type res:          :class:`ost.mol.ResidueView`
    :type env:          :class:`ost.mol.EntityView`
    :type normalize:    :class:`bool`

    :returns:           The requested energy
    :rtype:             :class:`float`

  .. method:: GetTotalEnergy(target, env, normalize)

    Gets the energy of **target** given the provided environment

    :param target:      Structure to be evaluated
    :param env:         The environment, the residues from **target** actually sees
    :normalize:         Whether the final summed energy should be normalized by
                        total number of observed interactions
    :type target:       :class:`ost.mol.EntityView`
    :type env:          :class:`ost.mol.EntityView`
    :type normalize:    :class:`bool`

    :returns:           The requested energy, NaN if no interaction could be observed
    :rtype:             :class:`float`

  .. method:: GetEnergies(target, env)

    Gets energies for every single residue in the **target** given the provided
    environment. 

    :param target:      Structure to be evaluated
    :param env:         The environment, the residues from **target** actually sees

    :type target:       :class:`ost.mol.EntityView`
    :type env:          :class:`ost.mol.EntityView`
    :type normalize:    :class:`bool`

    :returns:           The requested energies
    :rtype:             :class:`list` of :class:`float`

  .. method:: GetOptions()

    :returns:           A readonly option object
    :rtype:             :class:`CBPackingOpts`


.. _hbond-potential:

HBondPotential
--------------------------------------------------------------------------------


.. class:: HBondOpts

  Internal option class without init function. This class gets built when
  initializing the according statistic object and is automatically passed over
  when creating the potential. The class is considered to be readonly.

  .. attribute:: d_min

    Only HBonds with d > d_min are considered

  .. attribute:: d_max

    Only HBonds with d < d_max are considered

  .. attribute:: alpha_min

    Only HBonds with alpha > alpha_min are considered

  .. attribute:: alpha_max

    Only HBonds with alpha < alpha_max are considered

  .. attribute:: beta_min

    Only HBonds with beta > beta_min are considered

  .. attribute:: beta_max

    Only HBonds with beta < beta_max are considered


  .. attribute:: gamma_min

    Only HBonds with gamma > gamma_min are considered

  .. attribute:: gamma_max

    Only HBonds with gamma < gamma_max are considered

  .. attribute:: seq_sep

    Only HBonds from residues separated by this amount of residues
    are considered


.. class:: HBondStatistic(d_min, d_max, alpha_min, alpha_max, beta_min, 
                          beta_max, gamma_min, gamma_max, d_bins,
                          alpha_bins, beta_bins, gamma_bins, seq_sep)

  Builds histograms for each pair of the 3 states given the 
  input parameters that can be filled using structural data.

  :param d_min:         Minimal distance for HBond
  :param d_max:         Maximal distance for HBond
  :param alpha_min:     Minimal alpha angle in range [0,pi] for HBond
  :param alpha_max:     Maximal alpha angle in range [0,pi] for HBond
  :param beta_min:      Minimal beta angle in range [0,pi] for HBond
  :param beta_max:      Maximal beta angle in range [0,pi] for HBond
  :param gamma_min:     Minimal gamma dihedral angle in range [-pi,pi] for HBond
  :param gamma_max:     Maximal gamma dihedral angle in range [-pi,pi] for HBond
  :param d_bins:        Number of distance bins
  :param alpha_bins:    Number of alpha angle bins
  :param beta_bins:     Number of beta angle bins
  :param gamma_bins:    Number of gamma angle bins
  :param seq_sep:       Sequence separation parameter

  :type d_min:          :class:`float`
  :type d_max:          :class:`float`
  :type alpha_min:      :class:`float`
  :type alpha_max:      :class:`float`
  :type beta_min:       :class:`float`
  :type beta_max:       :class:`float`
  :type gamma_min:      :class:`float`
  :type gamma_max:      :class:`float`
  :type d_bins:         :class:`int`
  :type alpha_bins:     :class:`int`
  :type beta_bins:      :class:`int`
  :type gamma_bins:     :class:`int`
  :type seq_sep:        :class:`int`


  .. method:: Load(filename)

    Loads a :class:`HBondStatistic` from disk

    :param filename:    Path to file containing the statistic
    :type filename:     :class:`str`

    :returns:           The loaded statistic
    :rtype:             :class:`HBondStatistic`
    :raises:            :class:`RuntimeError` if **filename** does not exist

  .. method:: Save(filename)

    Saves the statistic down to disk

    :param filename:    Path of the file, the statistic should be saved in
    :type filename:     :class:`str`

  .. method:: Extract(target, [weight = 1.0])

    Adds the value of the **weight** parameter to the according bins in the
    internal histograms

    :param target:      The target structure from which the statistics are 
                        extracted
    :param weight:      The value that gets added to the particular histogram 
                        positions

    :type target:       :class:`ost.mol.EntityView`
    :type weight:       :class:`float`

  .. method:: GetTotalCount()

    Sums up all values in all histogram bins

    :returns:           The summed histogram value
    :rtype:             :class:`float`

  .. method:: GetCount(state)

    :param state:       The state... 0: not defined 1: helix 2: sheet

    :type state:        :class:`int`

    :returns:           The sum of all counts in the histogram for 
                        one particular state
    :rtype:             :class:`float`

  .. method:: GetCount(state, d_bin, alpha_bin, beta_bin, gamma_bin)

    :param state:       The state... 0: not defined 1: helix 2: sheet
    :param d_bin:       Distance bin
    :param alpha_bin:   Alpha angle bin
    :param beta_bin:    Beta angle bin
    :param gamma_bin:   Gamma dihedral angle bin

    :type state:        :class:`int`
    :type d_bin:        :class:`int`
    :type alpha_bin:    :class:`int`
    :type beta_bin:     :class:`int`
    :type gamma_bin:    :class:`int`

    :returns:           The sum of all counts in the histogram for
                        this particular bin combination

    :rtype:             :class:`float`

  .. method:: GetOptions()

    :returns:           A readonly option object
    :rtype:             :class:`HBondOpts`


.. class:: HBondPotential()

  Uses a :class:`HBondStatistic` to generate a potential. The object must be
  created with its static Create function.

  .. method:: Load(filename)

    Loads a :class:`HBondPotential` from disk

    :param filename:    Path to file containing the potential
    :type filename:     :class:`str`

    :returns:           The loaded potential
    :rtype:             :class:`HBondPotential`
    :raises:            :class:`RuntimeError` if **filename** does not exist

  .. method:: Save(filename)

    Saves a the potential down to disk

    :param filename:    Path of the file, the potential should be saved in
    :type filename:     :class:`str`  

  .. method:: Create(stat)

    Uses the passed statistic to generate probability density functions, that
    are turned into a potential using the formalism described by Sippl. Note, that
    there is no sigma parameter given, the statistics are expected to be saturated.

    :param stat:        The statistic from which you want to create the potential
    :param sigma:       The sigma parameter to handle sparse data

    :type stat:         :class:`HBondStatistic`

    :returns:           A newly created potential
    :rtype:             :class:`HBondPotential`    

  .. method:: GetEnergy(state, d, alpha, beta, gamma)

    Get an energy value for a particular set of parameters

    :param state:       The state... 0: not defined 1: helix 2: sheet
    :param d:           The distance
    :param alpha:       The alpha angle
    :param beta:        The beta angle
    :param gamma:       The gamma angle

    :type state:        :class:`int`
    :type d:            :class:`float`
    :type alpha:        :class:`float`
    :type beta:         :class:`float`
    :type gamma:        :class:`float`

    :returns:           The energy given the input parameters, NaN if one of the
                        parameters is out of the allowed range.
    :rtype:             :class:`float`

  .. method:: GetEnergy(res, env)

    Gets the energy of **res** given the provided environment.

    :param res:         Residue to be evaluated
    :param env:         The environment, the residue actually sees

    :type res:          :class:`ost.mol.ResidueView`
    :type env:          :class:`ost.mol.EntityView`

    :returns:           The requested energy, NaN if no interactions
                        can be observed
    :rtype:             :class:`float`

  .. method:: GetTotalEnergy(target, normalize)

    Gets the energy of **target**

    :param target:      Structure to be evaluated
    :normalize:         Whether the summed energy should be normalized by the
                        total number of residues in **target**

    :type target:       :class:`ost.mol.EntityView`
    :type normalize:    :class:`bool`

    :returns:           The requested energy, NaN if no interaction could be observed
    :rtype:             :class:`float`

  .. method:: GetEnergies(target)

    Gets energies for every single residue in the **target** 

    :param target:      Structure to be evaluated

    :type target:       :class:`ost.mol.EntityView`

    :returns:           The requested energies, elements can be NaN, if
                        no interactions are observed for a particular residue
    :rtype:             :class:`list` of :class:`float`

  .. method:: GetOptions()

    :returns:           A readonly option object
    :rtype:             :class:`HBondOpts`




.. [sippl1990] Sippl MJ (1990). Calculation of conformational ensembles from potentials of mean force. An approach to the knowledge-based prediction of local structures in globular proteins. J. Mol. Biol. 

.. [kortemme2003] Kortemme T, Morozov AV, Baker D (2003). An orientation-dependent hydrogen bonding potential improves prediction of specificity and structure for proteins and protein-protein complexes. J. Mol. Biol.
