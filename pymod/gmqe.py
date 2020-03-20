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

from ost import seq
from ost import mol
from ost import geom
from qmean import conf
from qmean import GMQEScoreCalculator
from qmean import PotentialContainer

class GMQE:

    def __init__(self, seqres, psipred, disco, 
                 profile=None, potentials = None):

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

        if str(disco.GetSeqres()) != str(self.seqres):
            raise RuntimeError('DisCoContainer is inconsistent with SEQRES')
        self.disco = disco

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


        self.score_calculator = GMQEScoreCalculator(self.potentials,
                                                    self.disco,
                                                    self.seqres,
                                                    self.psipred_pred,
                                                    self.psipred_cfi)


    def GetScores(self, aln, n_positions = None, ca_positions = None, 
                  c_positions = None, cb_positions = None, dssp_states = None):

        # expect exactly two sequences in aln
        if aln.GetCount() != 2:
            raise RuntimeError("Expect exactly two sequences in aln")

        # expect first sequence to match seqres 
        if str(aln.GetSequence(0).GetGaplessString()) != str(self.seqres):
            raise RuntimeError("Expect first seq in aln to match SEQRES")

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

        scores['seq_id'] = seq.alg.SequenceIdentity(aln)
        scores['seq_sim'] = seq.alg.SequenceSimilarity(aln, seq.alg.BLOSUM62, 
                                                       normalize=True)

        return scores

