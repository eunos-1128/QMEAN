# Copyright (c) 2013-2020, SIB - Swiss Institute of Bioinformatics and
# Biozentrum - University of Basel
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
# http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import re
from ost import seq
from ost import mol
from ost import geom
from qmean import DisCoContainer
from qmean import conf
from qmean import GMQEScoreCalculator
from qmean import PotentialContainer

class GMQE:

    def __init__(self, seqres, psipred, disco=None, 
                 profile=None, potentials = None, crf_file = None):

        if isinstance(seqres, str):
            self.seqres = seq.CreateSequence('seqres', seqres)
        else:
            self.seqres = seqres
        if '-' in str(self.seqres):
            raise RuntimeError("SEQRES must contain no gaps")

        if str(psipred.seq) != str(self.seqres):
            raise RuntimeError('PSIPREDHandler is inconsistent with SEQRES')
        # we're not directly using the QMEAN psipred scoring capabilities
        # but rather hack something together... lets store the actual psipred
        # prediction as a sequence that can then be fed into the 
        # GMQEScoreCalculator
        self.psipred_pred = seq.CreateSequence('psipred', ''.join(psipred.ss))
        self.psipred_cfi = [int(cfi) for cfi in psipred.conf]

        # make DisCoContainer optional... We still need to feed "some"
        # DisCoContainer into the GMQEScoreCalculator... If it's only an empty
        # container, the resulting score will be zero and useless. We therefore
        # keep track of that situation with the fake_disco variable.
        if disco is None:
            self.disco = DisCoContainer(self.seqres)
            self.fake_disco = True
        else:
            if str(disco.GetSeqres()) != str(self.seqres):
                raise RuntimeError('DisCoContainer is inconsistent with SEQRES')
            self.disco = disco
            self.fake_disco = False

        self._profile_with_pseudo_counts = None
        if profile:
            if str(profile.sequence) != str(self.seqres):
                raise RuntimeError('Profile is inconsistent with SEQRES')
            self.profile = profile
            self.profile_avg_entropy = profile.avg_entropy
        else:
            self.profile = None
            self.profile_avg_entropy = None

        if potentials:
            self.potentials = potentials
        else:
            s = conf.SwissmodelSettings()
            self.potentials = PotentialContainer.Load(s.global_potentials)

        self.crf_file = crf_file

        self.score_calculator = GMQEScoreCalculator(self.potentials,
                                                    self.disco,
                                                    self.seqres,
                                                    self.psipred_pred,
                                                    self.psipred_cfi)


    @property
    def profile_with_pseudo_counts(self):
        if self._profile_with_pseudo_counts is None:
            if self.profile is None:
                return None
            # enforce deep copy by using Extract function
            self._profile_with_pseudo_counts = self.profile.Extract(0, 
                                                         len(str(self.profile.sequence)))
            seq.alg.AddTransitionPseudoCounts(self._profile_with_pseudo_counts)
            if self.crf_file:
                crf_lib = seq.alg.ContextProfileDB.FromCRF(self.crf_file)
                seq.alg.AddAAPseudoCounts(self._profile_with_pseudo_counts, 
                                          crf_lib)
            else:
                seq.alg.AddAAPseudoCounts(self._profile_with_pseudo_counts)
            seq.alg.AddNullPseudoCounts(self._profile_with_pseudo_counts)

        return self._profile_with_pseudo_counts


    def GetScores(self, aln, seqres_aln, n_positions = None, ca_positions = None, 
                  c_positions = None, cb_positions = None, dssp_states = None,
                  profile_aln_score = None, tpl_profile = None):

        # expect exactly two sequences in aln
        if aln.GetCount() != 2:
            raise RuntimeError("Expect exactly two sequences in aln")

        # expect first sequence to match seqres 
        if str(aln.GetSequence(0).GetGaplessString()) != str(self.seqres):
            raise RuntimeError("Expect first seq in aln to match SEQRES")

        # expect exactly two sequences in seqres_aln
        if seqres_aln.GetCount() != 2:
            raise RuntimeError("Expect exactly two sequences in seqres_aln")

        # expect first sequence to match seqres 
        if str(seqres_aln.GetSequence(0).GetGaplessString()) != str(self.seqres):
            raise RuntimeError("Expect first seq in seqres_aln to match SEQRES")

        residue_numbers = list()
        if n_positions is None or ca_positions is None or c_positions is None or\
           cb_positions is None or dssp_states is None:
            # at least one of the positions is None, expect view to be attached
            # to sequence 1
            s = aln.GetSequence(1)
            if not s.HasAttachedView():
                raise RuntimeError("Not all positions and dssp states provided, "\
                                   "need to extract them from attached " +\
                                   "EntityView. " +\
                                   "But nothing is attached to second sequence")

            n_positions = geom.Vec3List()
            ca_positions = geom.Vec3List()
            c_positions = geom.Vec3List()
            cb_positions = geom.Vec3List()
            dssp_states = list()

            current_rnum = 0
            for col in aln:
                if col[0] != '-':
                    current_rnum += 1
                if col[0] != '-' and col[1] != '-':
                    r = col.GetResidue(1)
                    if r.IsValid(): 
                        n = r.FindAtom("N")
                        ca = r.FindAtom("CA")
                        c = r.FindAtom("C")
                        cb = r.FindAtom("CB")
                        if n.IsValid() and ca.IsValid() and c.IsValid():
                            dssp_states.append(str(r.GetSecStructure()))
                            residue_numbers.append(current_rnum)
                            n_positions.append(n.GetPos())
                            ca_positions.append(ca.GetPos())
                            c_positions.append(c.GetPos())
                            if cb.IsValid():
                                cb_positions.append(cb.GetPos())
                            else:
                                cb_pos = mol.alg.CBetaPosition(n.GetPos(),
                                                               ca.GetPos(),
                                                               c.GetPos())
                                cb_positions.append(cb_pos)

        else:
            current_rnum = 0
            for col in aln:
                if col[0] != '-':
                    current_rnum += 1
                if col[0] != '-' and col[1] != '-':
                    residue_numbers.append(current_rnum)

            # length consistency with all other data is checked in the
            # score_calculator
            if len(residue_numbers) != len(n_positions):
                raise RuntimeError("If directly providing structural data, " +\
                                   "length of all lists must exactly match " +\
                                   "number of aligned columns in aln!")


        # expect dssp states either as string or list of characters but
        # need ost.SequenceHandle for score_calculator
        d = seq.CreateSequence('dssp', ''.join([item for item in dssp_states]))

        scores = self.score_calculator.Eval(n_positions, ca_positions, 
                                            c_positions, cb_positions, 
                                            d, residue_numbers)

        # if disco was not given at initialization, "dist_const" in scores
        # is invalid and needs to be removed
        if self.fake_disco:
            del scores["dist_const"]

        scores['seq_id'] = seq.alg.SequenceIdentity(aln)
        scores['seq_sim'] = seq.alg.SequenceSimilarity(aln, seq.alg.BLOSUM62, 
                                                       normalize=True)
        scores['n_insertions'] = len(re.findall("-+", 
                                     str(aln.GetSequence(1)).strip('-')))
        scores['n_deletions'] = len(re.findall("-+", 
                                    str(aln.GetSequence(0)).strip('-')))
        scores['seqres_length'] = len(self.seqres.GetGaplessString())

        seqres_covered = 0
        seqres_tot = 0
        for col in seqres_aln:
            if col[0] != '-':
                seqres_tot += 1
                if col[1] != '-':
                    seqres_covered += 1
        scores["seqres_coverage"] = float(seqres_covered)/seqres_tot

        if profile_aln_score is not None:
            scores["profile_aln_score"] = profile_aln_score
            scores["avg_entropy"] = self.profile_avg_entropy
        elif self.profile is not None and tpl_profile is not None:
            tpl_profile_copy = tpl_profile.Extract(0, len(tpl_profile.sequence))
            seq.alg.AddTransitionPseudoCounts(tpl_profile_copy)
            seq.alg.AddAAPseudoCounts(tpl_profile_copy)
            seq.alg.AddNullPseudoCounts(tpl_profile_copy)
            scores["profile_aln_score"] = \
            seq.alg.HMMScore(self.profile_with_pseudo_counts, 
                             tpl_profile_copy, seqres_aln, 0, 1)
            scores["avg_entropy"] = self.profile_avg_entropy

        return scores
