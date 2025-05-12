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

Predicted Sequence Features
================================================================================

.. currentmodule:: qmean

QMEAN allows to evaluate the consistency of features predicted from sequence
with their actual outcome in the model. There is the consistency of the
predicted solvent accessibility by ACCPro and the consistency with the
secondary structure prediction by PSIPRED. Both Terms come with some kind
of Python wrapper to make the model specific handling easier. In the more
sophisticated case of the secondary structure agreement, there even 
exists a dedicated class for score calculation.



Secondary Structure Agreement
--------------------------------------------------------------------------------

.. literalinclude:: example_scripts/ss_agreement_example.py

.. class:: SSAgreement

  A class, that holds precalculated probabilities for the occurence of
  dssp and psipred states as well as the probability for them to occur
  together. It then allows to calculate a logodd score for a certain 
  combination of actual secondary structure in a model as calculated
  by dssp and the according PSIPRED prediction.

  .. method:: GetDSSPP(dssp_state)

    :param dssp_state:  The dssp state, one of ['H','E','C','G','B','S','T','I']

    :type dssp_state:   :class:`str`

    :returns:           The probability of **dssp_state**
    :rtype:             :class:`float`

    :raises:            :class:`RuntimeError` if **dssp_state** is invalid

  .. method:: GetPsiPredP(psipred_state, psipred_conf)

    :param psipred_state: The PSIPRED prediction, one of ['H','E','C']
    :param psipred_conf: The PSIPRED confidence, one of [0,1,2,3,4,5,6,7,8,9]

    :type psipred_state: :class:`str`
    :type psipred_conf: :class:`int`

    :returns:           The probability of observing **psipred_state** with the
                        according **psipred_confidence**
    :rtype:             :class:`float`

    :raises:            :class:`RuntimeError` if one of the parameters is
                        invalid

  .. method:: GetOverallP(dssp_state, psipred_state, psipred_conf)

    :param dssp_state:  The dssp state, one of ['H','E','C','G','B','S','T','I']
    :param psipred_state: The PSIPRED prediction, one of ['H','E','C']
    :param psipred_conf: The PSIPRED confidence, one of [0,1,2,3,4,5,6,7,8,9]

    :type dssp_state:   :class:`str`
    :type psipred_state: :class:`str`
    :type psipred_conf: :class:`int`

    :returns:           The probability of observing **psipred_state** with the
                        according **psipred_confidence** in combination with
                        **dssp_state**
    :rtype:             :class:`float`

    :raises:            :class:`RuntimeError` If one of the parameters is 
                        invalid

  .. method:: LogOddScore(dssp_state, psipred_state, psipred_conf)

    :param dssp_state:  The dssp state, one of ['H','E','C','G','B','S','T','I']
    :param psipred_state: The PSIPRED prediction, one of ['H','E','C']
    :param psipred_conf: The PSIPRED confidence, one of [0,1,2,3,4,5,6,7,8,9]

    :type dssp_state:   :class:`str`
    :type psipred_state: :class:`str`
    :type psipred_conf: :class:`int`

    :returns:           A logodd score => the logarithm of the fraction between
                        the probabilities of observing all the input parameters
                        in combination and observing them independently from
                        each other.

    :rtype:             :class:`float`

    :raises:            :class:`RuntimeError` If one of the parameters is 
                        invalid  


.. literalinclude:: example_scripts/psipred_handler_example.py


.. class:: PSIPREDHandler(data)

  Allows to relate a psipred prediction with structural data. The prediction
  is related to a SEQRES. The :class:`PSIPREDHandler` allow to relate this
  prediction with a structure given consistent sequence.

  :param data:          Input data for the PSIPRED prediction. Must be 
                        a dictionary with the mentioned elements
                        as single strings with keys "seq", "ss" and "conf"

  :type data:           :class:`dict`


  .. method:: GetPSIPREDSS([chain=None])

    :param chain:       A chain view with sequence consistent to the internal
                        SEQRES.

    :type chain:        :class:`ost.mol.ChainHandle` / :class:`ost.mol.ChainView`

    :returns:           If **chain** is None, the internal predicted 
                        secondary structure gets returned, otherwise the
                        **chain** gets aligned to the internal SEQRES and
                        cut accordingly, so the returned prediction matches
                        the residues in **chain**.
    :rtype:             :class:`str`

    :raises:            :class:`RuntimeError` If **chain** is not None and
                        its sequence cannot be aligned to the internal SEQRES


  .. method:: GetPSIPREDConf([chain=None])

    :param chain:       A chain view with sequence consistent to the internal
                        SEQRES.

    :type chain:        :class:`ost.mol.ChainHandle` / :class:`ost.mol.ChainView`

    :returns:           If **chain** is None, the internal PSIPRED confidence
                        gets returned, otherwise the
                        **chain** gets aligned to the internal SEQRES and
                        cut accordingly, so the returned confidences matches
                        the residues in **chain**.
    :rtype:             :class:`str`

    :raises:            :class:`RuntimeError` If **chain** is not None and
                        its sequence cannot be aligned to the internal SEQRES

  .. method:: GetSSAgreementFromChain(chain)

    :param chain:       A chain view with sequence consistent to the internal
                        SEQRES.

    :type chain:        :class:`ost.mol.ChainView` / :class:`ost.mol.ChainHandle`

    :returns:           SSAgreement values mapped onto the structure by
                        aligning the structure to the internal SEQRES.

    :rtype:             :class:`list` of :class:`float`

    :raises:            :class:`RuntimeError` if the chain cannot be aligned
                        to the internal SEQRES.           





Solvent Accessibility Agreement
--------------------------------------------------------------------------------


.. literalinclude:: example_scripts/accpro_handler_example.py

.. class:: ACCPROHandler(data)

  Allows to relate an accessibility prediction with structural data. The
  prediction is related to a SEQRES. The :class:`ACCPROHandler` allow to
  relate this prediction with a structure given consistent sequence.

  :param data:          Input data for the ACCPRO prediction. The 
                        required information is SEQRES and accessibility.
                        Must be a dictionary with the mentioned elements
                        as single strings with keys "seq", "acc".
                        The accessibility must be descrobed in the
                        letters 'b' (buried) and 'e' (exposed).

  :type data:           :class:`dict`


  .. method:: GetACCPROAccessibility([chain=None])

    :param chain:       A chain view with sequence consistent to the internal
                        SEQRES.

    :type chain:        :class:`ost.mol.ChainHandle` / :class:`ost.mol.ChainView`

    :returns:           If **chain** is None, the internal predicted 
                        accessibility gets returned, otherwise the
                        **chain** gets aligned to the internal SEQRES and
                        cut accordingly, so the returned prediction matches
                        the residues in **chain**.
    :rtype:             :class:`str`

    :raises:            :class:`RuntimeError` If **chain** is not None and
                        its sequence cannot be aligned to the internal SEQRES


  .. method:: GetACCAgreementFromChain(chain)

    :param chain:       A chain view with sequence consistent to the internal
                        SEQRES.

    :type chain:        :class:`ost.mol.ChainView` / :class:`ost.mol.ChainHandle`

    :returns:           ACCAgreement values mapped onto the structure by
                        aligning the structure to the internal SEQRES.

    :rtype:             :class:`list` of :class:`float`

    :raises:            :class:`RuntimeError` if the chain cannot be aligned
                        to the internal SEQRES.           



