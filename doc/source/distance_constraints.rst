Distance Constraints (DisCo)
================================================================================

.. currentmodule:: qmean

In general the computation of DisCo scores can be separated into two parts. First information from known homologous templates structures needs to be found, extracted and stored. In a second step the DisCo scores can be computed and finally they can be combined with the other terms of QMEAN to the QMEANDisCo score.



Getting Information from Homologous Templates
--------------------------------------------------------------------------------
The Distance Constraints rely on information about residue distances in template structures being homologous to the structure of interest. This information is put into :class:`DCData` objects. Templates can be found and clustered by functions provided in the class :class:`TemplateInformation`. The :class:`DCData` object can be created with the function :meth:`qmean.FillDCData`. 

.. class:: TemplateInformation(seqres, working_dir = '.')
  
  A class to manage information from homologous template structures

  :param seqres:        The SEQRES from the model of interest
  :param working_dir:   The working directory where the SWISS-MODEL project is
                        stored during computation 

  :type seqres:         :class:`str`
  :type working_dir:    :class:`str`


  .. method:: FindTemplates()

    Creates a SWISS-MODEL project and identifies homologous templates to the SEQRES. This function needs a running SWISS-MODEL installation and access to the SWISS-MODEL templates library.

    :returns:           Multiple sequence alignment of the SEQRES with homologous
                        templates and list with paths to the template structures in the SMTL. Paths are listed in the same order as the templates in the multiple sequence alignment.
    :rtype:             (:class:`ost.seq.AlignmentHandle`, :class:`list`)

  .. method:: ClusterTemplates(msaln)

    Clusters all the templates from the multiple sequence alignment by their pairwise sequence similarity 

    :param msaln:       Multiple sequence alignment of the SEQRES with homologous
                        templates. 
    :type msaln:        :class:`ost.seq.AlignmentHandle`


    :returns:           List containing indices of templates in **msaln** that
                        are in the same cluster (first template in **msaln** gets index 0)
    :rtype:             :class:`list` of :class:`list`


.. method:: FillDCData(msaln, cluster, filenames, subst,\
                       dist_cutoff = 15, exp_factor = 70 )

  Creates the DCData object

  :param msaln:         Multiple sequence alignment of the SEQRES with all the
                        identified homologous sequences
  :param cluster:       Indices in **msaln** of templates that are in the 
                        same cluster. Can be obtained from :meth:`qmean.TemplateInformation.ClusterTemplates`.
  :param filenames:     List with paths to the template structures. Structure files  
                        should have the format: *filename*.pdb.gz
  :param subst:         Substitution matrix to compute the sequence similarity 
                        (for example BLOSUM62)
  :param dist_cutoff:   Distance cutoff for pairwise residue distances in the
                        template structures. Only residues being closer than the specified distance cutoff will be considered  
  :param exp_factor:    Constant needed for the cluster weighting factors. When 
                        cluster weights are computed the average sequence similarity of the templates in the same cluster is multiplied with this constant.  
 

  :type msaln:          :class:`ost.seq.AlignmentHandle`
  :type cluster:        :class:`list` of :class:`list`
  :type filenames:      :class:`list`
  :type subst:          :class:`ost.seq.alg.SubstWeightMatrix`
  :type dist_cutoff:    :class:`int`
  :type exp_factor:     :class:`int`

  :returns:           Object with information from homologous templates
  :rtype:             :class:`qmean.DCData`



.. method:: LoadDCData(filename)

  Loads a DCData object from a file

  :param filename:      DCData object file

  :type filename:       :class:`str`

  :returns:           Object with information from homologous templates
  :rtype:             :class:`qmean.DCData`


.. method:: SaveDCData(data, filename)

  Saves a DCData object to a file

  :param data:          DCData object file
  :param filename:      Name of file to which DCData object should be saved

  :type data:           :class:`qmean.DCData`
  :type filename:       :class:`str`

  :returns:           Object with information from homologous templates
  :rtype:             :class:`qmean.DCData`


.. class:: DCData
  
  This class contains all the data necessary for the computation of the DisCo scores like residue distances in homologous template structures, information about template clusters and some constants.  

  .. attribute:: seq

    The SEQRES for which the DCData has been created

  .. attribute:: cluster_weights

    List of weights for the template clusters

  .. attribute:: cluster_seq_sim

    Average sequence similarities to the SEQRES of templates in a cluster 

  .. attribute:: cluster_seq_id

    Average sequence identities to the SEQRES of templates in a cluster 

  .. attribute:: dist_cutoff

    Distance cutoff for residue distances in the templates
    
  .. attribute:: exp_factor

    Constant for cluster weight computation 


.. literalinclude:: example_scripts/dcdata_computation_example.py


Computation of DisCo and QMEANDisCo Scores 
--------------------------------------------------------------------------------

Once the data from homologous templates of a certain SEQRES is prepared (:class:`qmean.DCData` object) the DisCo scores can be easily computed for all models of the SEQRES of interest. Furthermore all values needed for assembling the QMEAN and the DisCo scores can be determined.

.. method:: DCScore(aln, data , dist_cutoff=15)

  Computes the local DisCo score for each residue from the model of interest

  :param aln:           An alignment of the SEQRES with the sequence of the model 
                        to be assessed. Can be created with the ost function :meth:`ost.seq.alg.AlignToSEQRES`. The model structure has to be attached (with :meth:`ost.seq.AlignmentHandle.AttachView`)
 
  :param data:          Object with information from homologous templates
  :param dist_cutoff:   Distance cutoff for residue distances in  the model. Only
                        pairwise residue distances smaller than dist_cutoff will be incorporated in the score computation 

  :type aln:            :class:`ost.seq.AlignmentHandle`
  :type data:           :class:`DCData`
  :type dist_cutoff:    :class:`int`

  :returns:           List with DisCo scores for each residue of the model that 
                      is attached to **aln**
  :rtype:             :class:`list`


.. method:: DetermineFeatureValues(aln, data)

  :param aln:           An alignment of the SEQRES with the sequence of the model 
                        to be assessed. Can be created with the ost function :meth:`ost.seq.alg.AlignToSEQRES`. The model structure has to be attached (with :meth:`ost.seq.AlignmentHandle.AttachView`).
 
  :param data:          Object with information from homologous templates of the SEQRES 

  :type aln:            :class:`ost.seq.AlignmentHandle`
  :type data:           :class:`DCData`

  :returns:           List with feature values for each residue
  :rtype:             :class:`list` of :class:`list`

.. literalinclude:: example_scripts/qmeandisco_example.py