import argparse
import datetime
import json
import os
import hashlib
import math
import tempfile
import subprocess

import qmean
import ost
from ost import io
from ost import conop
from ost import mol
from ost import seq
from ost.bindings import hhblits3
from qmean.predicted_sequence_features import PSIPREDHandler
from qmean.predicted_sequence_features import ACCPROHandler
from qmean.predicted_sequence_features import AlignChainToSEQRES


class _MolckLogger(ost.LogSink):
    '''Append Molck output to results object.'''
    def __init__(self):
        ost.LogSink.__init__(self)
        self.full_message = ''

    def LogMessage(self, message, severity):
        self.full_message += message


def _get_seq_name(sequence):
    return hashlib.md5(str(sequence).encode()).hexdigest()


def _run_accpro(fasta_file, msa, workdir, seqlen):

    acc_in = os.path.join(workdir, 'accpro_in')
    acc_out = os.path.join(workdir, 'accpro_out')

    with open(acc_in, 'w') as fh:
        fh.write(f'{msa.GetCount()}\n')
        for s in msa.sequences:
            s_with_dot = str(s).replace('-', '.') 
            fh.write(f'{s_with_dot}\n')
    
    cmd = ['perl', '/qmean/predict_ss_sa.pl', fasta_file, acc_out, acc_in]
    cmd = ' '.join(cmd)
    job = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, 
                         stderr=subprocess.PIPE)
    
    sout = job.stdout.decode()
    serr = job.stderr.decode()
    if 'Segmentation fault' in sout:
        LogError('Running ACCPro failed for command "%s" with a segmentation fault' % command)
        return None
    if 'Segmentation fault' in serr:
        LogError('Running ACCPro failed for command "%s" with a segmentation fault' % command)
        return None
    if job.returncode != 0:
        LogError('Running ACCPro failed for command "%s":\n%s\n%s' %(' '.join(command),
                                                                            sout,
                                                                            serr))
        return None

    with open(acc_out, 'r') as fh:
        data = fh.readlines()
        if len(data) != 3:
            LogError('Running ACCPro failed, corrupt output')
            return None
        return data[2].strip().replace('-', 'b')


def _seqanno(s, workdir, uniclust30, do_disco=False):

    seqanno_workdir = os.path.join(workdir, s.GetName())
    hh = hhblits3.HHblits(s,  '/usr/local', working_dir=seqanno_workdir)
    a3m_file = hh.BuildQueryMSA(uniclust30)
    with open(a3m_file) as fh:
        a3m_content = hhblits3.ParseA3M(fh)

    hhm_file = hh.A3MToProfile(a3m_file)
    prof = hhblits3.ParseHHM(open(hhm_file))

    # data for PSIPREDHandler can be fetched from a3m_content
    data = dict()
    data['ss'] = a3m_content['ss_pred']
    data['conf'] = [int(c) for c in a3m_content['ss_conf']]
    data['seq'] = str(s)
    psipred_handler = PSIPREDHandler(data)

    # accpro must be called separately
    # accpro fails in rare circumstances we still want to continue in
    # this case (we do the same in SWISS-MODEL)
    acc = _run_accpro(hh.filename, prof['msa'], seqanno_workdir, len(s))
    if acc is None:
        accpro_handler = None
    else:
        data['acc'] = acc
        accpro_handler = ACCPROHandler(data)

    if do_disco:
        dc = None # do fancy stuff
    else:
        dc = None

    return psipred_handler, accpro_handler, dc



class ModelScorer:
    def __init__(self, model_path, seqres):
        self.model_path = model_path
        self.seqres = seqres

        # sets the following attributes:
        # - model: Model loaded from *model_path*
        self._load_model()

        # sets the following attributes:
        # - processed_model: model with Molck applied (see code) and potentially
        #                    added chain name.
        #                    Residues are renumbered according to SEQRES
        # - chain_name_mapping: Dictionary mapping chain names in model to 
        #                       processed_model
        # - seqres_list: List of length 
        #                len(self.processed_model.Select('peptide=true').chains)
        #                Contains for each chain the SEQRES sequence. 
        # - atomseq_list: List of length 
        #                 len(self.processed_model.Select('peptide=true').chains)
        #                 Contains for each chain the ATOMSEQ sequence. 
        self._process_model()

        # the following members remain empty until you call score()
        self.local_scores = None
        self.global_scores = None


    def to_json(self):
        out_dict = dict()
        out_dict['chains'] = dict()
        for ch, atomseq, seqres in zip(self.processed_model.chains, 
                                       self.atomseq_list, self.seqres_list):
            out_dict['chains'][ch.GetName()] = {'atomseq': str(atomseq),
                                                'seqres': str(seqres)}
        out_dict['original_name'] = self.model_path
        out_dict['preprocessing'] = self.preprocessing_log
        out_dict['scores'] = {'local_scores': self.local_scores,
                              'global_scores': self.global_scores}
        return out_dict


    def score(self, psipred_handler, accpro_handler, disco_container,
              scoring_function):

        psipred_handler_list = list()
        accpro_handler_list = list()
        disco_container_list = None

        for s in self.seqres_list:
            psipred_handler_list.append(psipred_handler[s.GetName()])
            accpro_handler_list.append(accpro_handler[s.GetName()])

        scorer = qmean.QMEANScorer(self.processed_model.Select('peptide=True'), 
                                   psipred=psipred_handler_list, 
                                   accpro=accpro_handler_list, 
                                   dc=disco_container_list, 
                                   use_nn=scoring_function=="QMEANDisCo")

        # QMEAN gives NaN if something is not defined. JSON prefers None
        def NaN_To_None(val):
            if math.isnan(val):
                return None
            return val

        global_scores = dict()
        global_scores["qmean4_norm_score"] = NaN_To_None(scorer.qmean4_score)
        global_scores["qmean4_z_score"] = NaN_To_None(scorer.qmean4_z_score)
        try:
            global_scores["qmean6_norm_score"] = NaN_To_None(scorer.qmean6_score)
            global_scores["qmean6_z_score"] = NaN_To_None(scorer.qmean6_z_score)
        except Exception:  # pylint: disable=broad-except
            global_scores["qmean6_norm_score"] = None
            global_scores["qmean6_z_score"] = None
        scores = scorer.qmean6_components
        global_scores["interaction_norm_score"] = NaN_To_None(scores["interaction"])
        global_scores["interaction_z_score"] = NaN_To_None(scores["interaction_z_score"])
        global_scores["cbeta_norm_score"] = NaN_To_None(scores["cbeta"])
        global_scores["cbeta_z_score"] = NaN_To_None(scores["cbeta_z_score"])
        global_scores["packing_norm_score"] = NaN_To_None(scores["packing"])
        global_scores["packing_z_score"] = NaN_To_None(scores["packing_z_score"])
        global_scores["torsion_norm_score"] = NaN_To_None(scores["torsion"])
        global_scores["torsion_z_score"] = NaN_To_None(scores["torsion_z_score"])
        global_scores["ss_agreement_norm_score"] = NaN_To_None(scores["ss_agreement"])
        global_scores["ss_agreement_z_score"] = NaN_To_None(scores["ss_agreement_z_score"])
        global_scores["acc_agreement_norm_score"] = NaN_To_None(scores["acc_agreement"])
        global_scores["acc_agreement_z_score"] = NaN_To_None(scores["acc_agreement_z_score"])
        global_scores["avg_local_score"] = NaN_To_None(scorer.avg_local_score)
        global_scores["avg_local_score_error"] = NaN_To_None(scorer.avg_local_score_error)

        local_scores = dict()
        for chn, seq in zip(scorer.model.chains, self.seqres_list):
            score_list = list([None] * len(seq))
            per_chain_scores = scorer.local_scores[chn.GetName()]
            for rnum, score in per_chain_scores.items():
                if rnum < 1 or rnum > len(seq):
                    raise RuntimeError(
                        "Observed ResNum not matching provided SEQRES!"
                    )
                score_list[rnum - 1] = NaN_To_None(score)
            local_scores[chn.GetName()] = score_list

        self.local_scores = local_scores
        self.global_scores = global_scores


    def _load_model(self):
        """Read OST entity either from mmCIF or PDB."""

        # Determine file format from suffix.
        ext = self.model_path.split(".")
        if ext[-1] == "gz":
            ext = ext[:-1]
        if len(ext) <= 1:
            raise RuntimeError(f'Could not determine format of file {self.model_path}.')
        sformat = ext[-1].lower()

        # Load the structure
        if sformat in ["mmcif", "cif"]:
            entity = io.LoadMMCIF(self.model_path)
            if len(entity.residues)==0:
                raise Exception(f'No residues found in model file: {self.model_path}')
        elif sformat in ['pdb']:
            entity = io.LoadPDB(self.model_path)
            if len(entity.residues)==0:
                raise Exception(f'No residues found in model file: {self.model_path}')
        else:
            raise RuntimeError(f'Unknow/ unsuppoted file extension found for file {self.model_path}.')
        self.model = entity


    def _process_model(self):

        # perform processing on deep copy of input model
        self.processed_model = self.model.Copy()

        #########################################################
        # Add chain names if not set (happens with CASP models) #
        #########################################################
        stupid_names = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_-"
        ed = self.processed_model.handle.EditXCS()
        self.chain_name_mapping = dict()
        for ch in self.processed_model.chains:
            if len(ch.GetName().strip()) < 1:
                # iterate over possible chain names and assign first unused name
                for new_cname in stupid_names:
                    try:
                        old_cname = ch.GetName()
                        ed.RenameChain(ch, new_cname)
                        self.chain_name_mapping[old_cname] = new_cname
                        break
                    except:
                        continue
            else:
                # chain name is valid, nothing to do for this chain
                self.chain_name_mapping[ch.GetName()] = ch.GetName()

        ############################
        # Preprocessing with Molck #
        ############################
        pre_atom_count = len(self.processed_model.atoms)
        pre_res_count = len(self.processed_model.residues)
        ost.PushVerbosityLevel(ost.LogLevel.Info)
        molcklogger = _MolckLogger()
        ost.PushLogSink(molcklogger)
        molck_settings = mol.alg.MolckSettings(rm_unk_atoms=True, 
                                               rm_non_std=True,
                                               rm_hyd_atoms=True, 
                                               rm_oxt_atoms=False,
                                               rm_zero_occ_atoms=True, 
                                               map_nonstd_res=True,
                                               assign_elem=True)

        mol.alg.Molck(self.processed_model, conop.GetDefaultLib(), 
                      molck_settings, prune=True)
        ost.PopLogSink()
        ost.PopVerbosityLevel()
        self.preprocessing_log = dict()
        self.preprocessing_log['molck_log'] = molcklogger.full_message
        self.preprocessing_log['atomsRemoved'] = pre_atom_count - len(self.processed_model.atoms)
        self.preprocessing_log['residuesRemoved'] = pre_res_count - len(self.processed_model.residues)
        self.peptide_processed_model = self.processed_model.Select('peptide=true')

        #########################
        # Etract SEQRES/ATOMSEQ #
        #########################
        alignments = list()
        if self.seqres:
            # SEQRES is provided by user, requires mapping to model chains

            # option 1: all chains align to this single SEQRES (monomer or homo-oligomer)
            if len(self.seqres) == 1:
                for ch in self.peptide_processed_model:
                    try:
                        aln = AlignChainToSEQRES(ch, self.seqres[0])
                        alignments.append(aln)
                    except:
                        raise RuntimeError(f'Failed to align chain {ch.GetName()} of {self.model_path} to provided SEQRES.')

            # option2: Map chains using names (whatever-mer, as long as
            #          as I find a SEQRES for each chain based on name matching)
            elif len(self.seqres) > 1:
                for ch in self.peptide_processed_model:
                    ch_seqres = None
                    for s in seqres:
                        if ch.GetName() == s.GetName():
                            ch_seqres = s
                            break
                    if ch_seqres is None:
                        raise RuntimeError(f'Failed to find SEQRES for chain of name {ch.GetName()} in provided SEQRES list.')
                    try:
                        aln = AlignChainToSEQRES(ch, ch_seqres)
                        alignments.append(aln)
                    except:
                        raise RuntimeError(f'Failed to align chain {ch.GetName()} of {self.model_path} to provided SEQRES.')
        else:
            # No SEQRES provided, extract SEQRES from protein => SEQRES==ATOMSEQ
            for ch in self.peptide_processed_model.chains:
                s = ''.join([r.one_letter_code for r in ch.residues])
                aln = seq.CreateAlignment()
                aln.AddSequence(seq.CreateSequence('SEQRES', s))
                aln.AddSequence(seq.CreateSequence('ATOMSEQ', s))
                alignments.append(aln)
        self.seqres_list = list()
        self.atomseq_list = list()
        for aln in alignments:
            s = aln.GetSequence(0).GetGaplessString()
            a = aln.GetSequence(1).GetGaplessString()
            self.seqres_list.append(seq.CreateSequence(_get_seq_name(s), s))
            self.atomseq_list.append(seq.CreateSequence(_get_seq_name(a), a))


        ################################
        # Renumber according to SEQRES #
        ################################
        ed = self.peptide_processed_model.handle.EditXCS()
        for aln, ch in zip(alignments, self.peptide_processed_model.chains):
            view = self.peptide_processed_model.CreateEmptyView()
            view.AddChain(ch, mol.INCLUDE_ALL)
            aln.AttachView(1, view)
            for i in range(aln.GetLength()):
                res = aln.GetResidue(1, i)
                if res.IsValid():
                    ed.SetResidueNumber(res.handle, ost.mol.ResNum(i + 1))


class ModelScorerContainer:
    def __init__(self, method, model_paths, seqres_path, workdir, uniclust30):
        self.created = datetime.datetime.now().isoformat(timespec='seconds')
        self.method = method
        self.qmean_version = qmean.qmean_version
        self.smtl_version = None
        self._seqres_list = None

        # load SEQRES and setup models
        seqres = self._load_seqres(seqres_path)
        self.models = list()
        for mp in model_paths:
            self.models.append(ModelScorer(mp, seqres))

        # set SEQRES specific data (dicts with key: seqres name, value: data):
        # - self.psipred_handler
        # - self.accpro_handler
        # - self.disco_container
        self._seqanno(workdir, uniclust30)

        # do scoring
        self.score()


    def score(self):
        for m in self.models:
            m.score(self.psipred_handler, self.accpro_handler, 
                    self.disco_container, self.method)


    def to_json(self, out_path):
        out_dict = dict()
        out_dict['created'] = self.created
        out_dict['method'] = self.method
        out_dict['qmean_version'] = self.qmean_version
        out_dict['smtl_version'] = self.smtl_version
        out_dict['models'] = dict()
        for m_idx, m in enumerate(self.models):
            identifier = str(m_idx)
            identifier = (3-len(identifier))*'0' + identifier
            out_dict['models'][identifier] = m.to_json()
        return out_dict


    @property
    def seqres_list(self):
        if not self._seqres_list:
            self._seqres_list = list()
            added_sequences = set()
            for m in self.models:
                for s in m.seqres_list:
                    if s.GetName() not in added_sequences:
                        self._seqres_list.append(s)
                        added_sequences.add(s.GetName())
        return self._seqres_list


    def _load_seqres(self, seqres_path):
        seqres = None
        if seqres_path:
            try:
                seqres = io.LoadSequenceList(seqres_path)
            except:
                raise RuntimeError(f'failed to load seqres file: {seqres_path}')
        return seqres


    def _seqanno(self, workdir, uniclust30):
        self.psipred_handler = dict()
        self.accpro_handler = dict()
        self.disco_container = dict()
        for s in self.seqres_list:
            p,a,d = _seqanno(s, workdir, uniclust30, 
                             do_disco = self.method == "QMEANDisCo")
            self.psipred_handler[s.GetName()] = p
            self.accpro_handler[s.GetName()] = a
            self.disco_container[s.GetName()] = d



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--method', dest='method', choices=['QMEAN', 
                        'QMEANDisCo', 'QMEANBrane'], default='QMEANDisCo')
    parser.add_argument('--models', required=True, dest='models', nargs='+')
    parser.add_argument('--out', dest='out', default='out.json')
    parser.add_argument('--seqres', dest='seqres', default=None)
    parser.add_argument('--workdir', dest='workdir', required=True)
    parser.add_argument('--complib', dest='complib', required=True)
    parser.add_argument('--uniclust30', dest='uniclust30', required=True)
    args = parser.parse_args()


    ##########################
    # Setup and input checks #
    ##########################
    for model_path in args.models:
        if not os.path.exists(model_path):
            raise RuntimeError(f'specified path {model_path} does not exist')

    if args.seqres:
        if not os.path.exists(args.seqres):
            raise RuntimeError(f'specified path {args.seqres} does not exist')

    if not os.path.exists(args.complib):
        raise RuntimeError(f'specified path {args.complib} does not exist')

    expected_uniclust30_suffixes = ['_a3m.ffdata', '_a3m.ffindex',
                                    '_hhm.ffdata', '_hhm.ffindex',
                                    '_cs219.ffdata', '_cs219.ffindex']
    for suffix in expected_uniclust30_suffixes:
        full = args.uniclust30 + suffix
        if not os.path.exists(full):
            raise RuntimeError(f'Expect {full} to be present in uniclust30')

    # load and set user defined compound library, don't assume it to be available
    # from the container
    if not os.path.exists(args.complib):
        raise RuntimeError(f'specified path {args.complib} does not exist')
    complib = conop.CompoundLib.Load(args.complib)
    conop.SetDefaultLib(complib)


    ######################################
    # Load/process models and do scoring #
    ######################################
    data = ModelScorerContainer(args.method, args.models, args.seqres, 
                                args.workdir, args.uniclust30)


    ###############
    # Dump output #
    ###############
    with open(args.out, 'w') as fh:
        json.dump(data.to_json(args.out), fh)


if __name__ == "__main__":
    main()

