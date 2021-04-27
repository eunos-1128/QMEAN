#!/usr/local/bin/ost

import argparse
import datetime
import json
import os
import hashlib
import math
import tempfile
import subprocess
import sys
import shutil

import ost
from ost import io
from ost import conop
from ost import mol
from ost import seq
from ost.bindings import hhblits3

import qmean
import qmean.conf as qmean_config
from qmean.predicted_sequence_features import PSIPREDHandler
from qmean.predicted_sequence_features import ACCPROHandler
from qmean.predicted_sequence_features import AlignChainToSEQRES
from qmean import DisCoContainer
from qmean import ExtractTemplateDataDisCo
from qmean import DisCoDataContainers
from qmean import mqa_result_membrane


class _Logger(ost.LogSink):
    """Append Molck output to results object."""

    def __init__(self):
        ost.LogSink.__init__(self)
        self.full_message = ""
        self.error_message = ""
        self.info_message = ""

    def LogMessage(self, message, severity):
        self.full_message += message
        # deal with errors
        if severity == ost.LogLevel.Error:
            self.error_message += message
        # deal with info
        if severity == ost.LogLevel.Info:
            self.info_message += message


class _ChainClusterIndex:
    def __init__(self, index_path, dates_path):
        # both are lazy loaded, in particular the dates might never be used...
        self.index_path = index_path
        self.dates_path = dates_path
        self._index = None
        self._dates = None

    def get_clusters(self, seq_hash, datefilter=None):

        clusters = list()
        if seq_hash in self.index:
            clusters = self.index[seq_hash]

        if datefilter is None:
            return clusters

        filtered_clusters = list()
        thresh_date = datetime.datetime.strptime(datefilter, "%Y-%m-%d")
        for cluster in clusters:
            filtered_cluster = list()
            for smt_id in cluster:
                rel_date = datetime.datetime.strptime(
                    self.get_date(smt_id), "%Y-%m-%d"
                )
                if rel_date < thresh_date:
                    filtered_cluster.append(smt_id)
            if len(filtered_cluster) > 0:
                filtered_clusters.append(filtered_cluster)
        return filtered_clusters

    def get_date(self, smt_id):
        pdb_id = smt_id.split(".")[0].upper()
        if pdb_id not in self.dates:
            raise RuntimeError(f"Could not find {pdb_id} in dates file")
        return self.dates[pdb_id]

    @property
    def index(self):
        if self._index is None:
            self._index = dict()
            with open(self.index_path, "r") as fh:
                stuff = fh.readlines()
            for line in stuff:
                split_line = line.split()
                clusters = list()
                for item in split_line[1:]:
                    clusters.append(item.split(","))
                self._index[split_line[0]] = clusters
        return self._index

    @property
    def dates(self):
        if self._dates is None:
            with open(self.dates_path, "r") as fh:
                stuff = fh.readlines()
            header = stuff[0].split(",")
            if "pdb_id" not in header or "release_date" not in header:
                raise RuntimeError(f"ill formatted file: {self.dates_path}")
            pdb_id_idx = header.index("pdb_id")
            release_date_idx = header.index("release_date")
            for line in stuff[1:]:
                split_line = line.split(",")
                pdb_id = split_line[pdb_id_idx].upper()
                date = split_line[release_date_idx]
                self._dates[pdb_id] = date
        return self._dates


def _get_seq_name(sequence):
    return hashlib.md5(str(sequence).encode()).hexdigest()


def _run_accpro(fasta_file, msa, workdir, seqlen):

    acc_in = os.path.join(workdir, "accpro_in")
    acc_out = os.path.join(workdir, "accpro_out")

    with open(acc_in, "w") as fh:
        fh.write(f"{msa.GetCount()}\n")
        for s in msa.sequences:
            s_with_dot = str(s).replace("-", ".")
            fh.write(f"{s_with_dot}\n")

    cmd = ["perl", "/qmean/predict_ss_sa.pl", fasta_file, acc_out, acc_in]
    job = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    sout = job.stdout.decode()
    serr = job.stderr.decode()
    if "Segmentation fault" in sout:
        LogError(
            'Running ACCPro failed for command "%s" with a segmentation fault'
            % command
        )
        return None
    if "Segmentation fault" in serr:
        LogError(
            'Running ACCPro failed for command "%s" with a segmentation fault'
            % command
        )
        return None
    if job.returncode != 0:
        LogError(
            'Running ACCPro failed for command "%s":\n%s\n%s'
            % (" ".join(command), sout, serr)
        )
        return None

    with open(acc_out, "r") as fh:
        data = fh.readlines()
        if len(data) != 3:
            LogError("Running ACCPro failed, corrupt output")
            return None
        return data[2].strip().replace("-", "b")


def _get_dc(
    workdir, target_seq, hh, a3m_file, smtldir, datefilter=None, store_dc=True
):

    db = os.path.join(smtldir, "smtl_uniq")

    # hhblits prints stuff in stderr, let's setup yet another logger and swallow
    # everything. Check later if something went wrong...
    logger = _Logger()
    ost.PushLogSink(logger)
    hhr_file = hh.Search(a3m_file, db, prefix="hhm100", options={"premerge": 0})
    ost.PopLogSink()

    if hhr_file is None:
        LogError("HHblits search did not produce any results, aborting")
        sys.exit(-1)

    index_file = os.path.join(smtldir, "CHAINCLUSTERINDEX")
    dates_file = os.path.join(smtldir, "dates.csv")
    chain_cluster_index = _ChainClusterIndex(index_file, dates_file)

    binary_tpl = DisCoDataContainers(
        os.path.join(smtldir, "indexer.dat"),
        os.path.join(smtldir, "seqres_data.dat"),
        os.path.join(smtldir, "atomseq_data.dat"),
        os.path.join(smtldir, "ca_pos_data.dat"),
    )

    dc = DisCoContainer(target_seq)
    header, hits = hhblits3.ParseHHblitsOutput(open(hhr_file))
    target_len = len(target_seq)

    for hit in hits:
        start = hit.aln.sequences[0].offset
        end = start + len(hit.aln.sequences[0].gapless_string)
        clusters = chain_cluster_index.get_clusters(hit.hit_id, datefilter)
        for cluster in clusters:
            tpl = cluster[0]
            aln = ost.seq.CreateAlignment(
                target_seq.Copy(), ost.seq.CreateSequence(tpl, "-" * target_len)
            )
            aln[start:end].Replace(hit.aln[:])
            aln.SetSequenceOffset(1, hit.aln.sequences[1].offset)
            pdb_id, assembly_id, chain_name = tpl.split(".")

            # use try except block, if certain template cannot be
            # found in binary SMTL...
            try:
                binary_tpl_assembly = ".".join([pdb_id, assembly_id])
                tpl_data = ExtractTemplateDataDisCo(
                    binary_tpl_assembly, chain_name, aln, binary_tpl
                )

            except Exception as e:
                ost.LogError(f"Cannot extract data from template {tpl} ({e})")
                continue
            # THE ALIGNMENT SHOULD BE THE ATOMSEQ ALN AS IN A SWISS-MODEL PROJECT
            dc.AddData(aln, tpl_data.ca_positions, tpl_data.residue_numbers)

    # just set everything to default values...
    disco_dist_thresh = 15.0
    disco_gamma = 70.0
    disco_seqsim_thresh = 0.5
    disco_bin_size = 0.5
    dc.CalculateConstraints(
        disco_dist_thresh, disco_gamma, disco_seqsim_thresh, disco_bin_size
    )

    if store_dc:
        out_path = os.path.join(workdir, "qmean_dc.dat")
        dc.Save(out_path)

    return dc


def _seqanno(target_seq, workdir, uniclust30, do_disco, smtldir, datefilter):

    md5_hash = target_seq.GetName()
    seqanno_workdir = os.path.join(workdir, md5_hash)
    hh = hhblits3.HHblits(target_seq, "/usr/local", working_dir=seqanno_workdir)

    # in case of user specified profiles or when a workdir has been defined
    # which was not cleaned up, there already is an a3m file
    a3m = os.path.join(seqanno_workdir, hhblits3.HHblits.OUTPUT_PREFIX + ".a3m")
    a3m_content = None
    if os.path.exists(a3m):
        ost.LogInfo(f"Use precomputed a3m file (SEQRES md5 hash: {md5_hash})")
        with open(a3m) as fh:
            a3m_content = hhblits3.ParseA3M(fh)
    else:
        # hhblits prints stuff in stderr, let's setup yet another logger and
        # swallow everything. Check later if something went wrong...
        logger = _Logger()
        ost.PushLogSink(logger)
        a3m = hh.BuildQueryMSA(uniclust30)
        with open(a3m) as fh:
            a3m_content = hhblits3.ParseA3M(fh)
        ost.PopLogSink()

    # hhblits prints stuff in stderr, let's setup yet another logger and swallow
    # everything. Check later if something went wrong...
    logger = _Logger()
    ost.PushLogSink(logger)
    hhm = hh.A3MToProfile(a3m)
    ost.PopLogSink()
    prof = hhblits3.ParseHHM(open(hhm))

    # data for PSIPREDHandler can be fetched from a3m_content
    data = dict()
    data["ss"] = a3m_content["ss_pred"]
    data["conf"] = [int(c) for c in a3m_content["ss_conf"]]
    data["seq"] = str(target_seq)
    psipred_handler = PSIPREDHandler(data)

    # accpro must be called separately
    # accpro fails in rare circumstances we still want to continue in
    # this case (we do the same in SWISS-MODEL)
    acc = _run_accpro(
        hh.filename, prof["msa"], seqanno_workdir, len(target_seq)
    )
    if acc is None:
        accpro_handler = None
    else:
        data["acc"] = acc
        accpro_handler = ACCPROHandler(data)

    if do_disco:
        dc = _get_dc(seqanno_workdir, target_seq, hh, a3m, smtldir, datefilter)
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
        out_dict["chains"] = dict()
        for ch, atomseq, seqres in zip(
            self.processed_model.chains, self.atomseq_list, self.seqres_list
        ):
            out_dict["chains"][ch.GetName()] = {
                "atomseq": str(atomseq),
                "seqres": str(seqres),
            }
        out_dict["original_name"] = self.model_path
        out_dict["preprocessing"] = self.preprocessing_log
        out_dict["scores"] = {
            "local_scores": self.local_scores,
            "global_scores": self.global_scores,
        }
        return out_dict

    def score(
        self, psipred_handler, accpro_handler, disco_container, scoring_function
    ):

        psipred_handler_list = list()
        accpro_handler_list = list()
        disco_container_list = None

        for s in self.seqres_list:
            psipred_handler_list.append(psipred_handler[s.GetName()])
            accpro_handler_list.append(accpro_handler[s.GetName()])

        if scoring_function == "QMEANDisCo":
            disco_container_list = list()
            for s in self.seqres_list:
                disco_container_list.append(disco_container[s.GetName()])

        scorer = qmean.QMEANScorer(
            self.processed_model.Select("peptide=True"),
            psipred=psipred_handler_list,
            accpro=accpro_handler_list,
            dc=disco_container_list,
            use_nn=scoring_function == "QMEANDisCo",
        )

        # QMEAN gives NaN if something is not defined. JSON prefers None
        def format_float(val, digits=3):
            if math.isnan(val):
                return None
            return round(val, digits)

        global_scores = dict()
        global_scores["qmean4_norm_score"] = format_float(scorer.qmean4_score)
        global_scores["qmean4_z_score"] = format_float(scorer.qmean4_z_score)
        try:
            global_scores["qmean6_norm_score"] = format_float(
                scorer.qmean6_score
            )
            global_scores["qmean6_z_score"] = format_float(
                scorer.qmean6_z_score
            )
        except Exception:  # pylint: disable=broad-except
            global_scores["qmean6_norm_score"] = None
            global_scores["qmean6_z_score"] = None
        scores = scorer.qmean6_components
        global_scores["interaction_norm_score"] = format_float(
            scores["interaction"]
        )
        global_scores["interaction_z_score"] = format_float(
            scores["interaction_z_score"]
        )
        global_scores["cbeta_norm_score"] = format_float(scores["cbeta"])
        global_scores["cbeta_z_score"] = format_float(scores["cbeta_z_score"])
        global_scores["packing_norm_score"] = format_float(scores["packing"])
        global_scores["packing_z_score"] = format_float(
            scores["packing_z_score"]
        )
        global_scores["torsion_norm_score"] = format_float(scores["torsion"])
        global_scores["torsion_z_score"] = format_float(
            scores["torsion_z_score"]
        )
        global_scores["ss_agreement_norm_score"] = format_float(
            scores["ss_agreement"]
        )
        global_scores["ss_agreement_z_score"] = format_float(
            scores["ss_agreement_z_score"]
        )
        global_scores["acc_agreement_norm_score"] = format_float(
            scores["acc_agreement"]
        )
        global_scores["acc_agreement_z_score"] = format_float(
            scores["acc_agreement_z_score"]
        )
        global_scores["avg_local_score"] = format_float(scorer.avg_local_score)
        global_scores["avg_local_score_error"] = format_float(
            scorer.avg_local_score_error
        )

        local_scores = dict()

        if scoring_function in ["QMEANDisCo", "QMEAN"]:
            for chn, seq in zip(scorer.model.chains, self.seqres_list):
                score_list = list([None] * len(seq))
                per_chain_scores = scorer.local_scores[chn.GetName()]
                for rnum, score in per_chain_scores.items():
                    if rnum < 1 or rnum > len(seq):
                        raise RuntimeError(
                            "Observed ResNum not matching provided SEQRES!"
                        )
                    score_list[rnum - 1] = format_float(score)
                local_scores[chn.GetName()] = score_list
        elif scoring_function == "QMEANBrane":
            # the global scores are the same as QMEAN but the local ones change
            settings = qmean_config.MembraneSettings()
            peptide_sel = self.processed_model.Select("peptide=True")
            res = mqa_result_membrane.LocalMembraneResult.Create(
                peptide_sel,
                settings,
                False,
                psipred=psipred_handler_list,
                accpro=accpro_handler_list,
            )

            for ch, s in zip(peptide_sel.chains, self.seqres_list):
                score_list = list([None] * len(s))
                local_scores[ch.GetName()] = score_list

            chain_name_idx = res.score_table.GetColIndex("chain")
            rnum_idx = res.score_table.GetColIndex("rnum")
            score_idx = res.score_table.GetColIndex("QMEAN")

            for r in res.score_table.rows:
                chain_name = r[chain_name_idx]
                rnum = r[rnum_idx]
                score = r[score_idx]
                local_scores[chain_name][rnum - 1] = format_float(score)
        else:
            raise RuntimeError(f"Unknown scoring function {scoring_function}")

        self.local_scores = local_scores
        self.global_scores = global_scores

    def _load_model(self):
        """Read OST entity either from mmCIF or PDB."""

        # Determine file format from suffix.
        ext = self.model_path.split(".")
        if ext[-1] == "gz":
            ext = ext[:-1]
        if len(ext) <= 1:
            raise RuntimeError(
                f"Could not determine format of file {self.model_path}."
            )
        sformat = ext[-1].lower()

        # increase loglevel, as we would pollute the info log with weird stuff
        ost.PushVerbosityLevel(ost.LogLevel.Error)
        # Load the structure
        if sformat in ["mmcif", "cif"]:
            entity = io.LoadMMCIF(self.model_path)
            if len(entity.residues) == 0:
                raise Exception(
                    f"No residues found in model file: {self.model_path}"
                )
        elif sformat in ["pdb"]:
            entity = io.LoadPDB(self.model_path)
            if len(entity.residues) == 0:
                raise Exception(
                    f"No residues found in model file: {self.model_path}"
                )
        else:
            raise RuntimeError(
                f"Unknown/ unsupported file extension found for file {self.model_path}."
            )

        # restore old loglevel
        ost.PopVerbosityLevel()
        self.model = entity

    def _process_model(self):

        self.preprocessing_log = dict()

        # perform processing on deep copy of input model
        self.processed_model = self.model.Copy()

        ############################
        # Preprocessing with Molck #
        ############################
        pre_atom_count = len(self.processed_model.atoms)
        pre_res_count = len(self.processed_model.residues)
        pre_ch_count = len(self.processed_model.chains)
        ost.PushVerbosityLevel(ost.LogLevel.Info)
        molcklogger = _Logger()
        ost.PushLogSink(molcklogger)
        molck_settings = mol.alg.MolckSettings(
            rm_unk_atoms=True,
            rm_non_std=True,
            rm_hyd_atoms=True,
            rm_oxt_atoms=False,
            rm_zero_occ_atoms=True,
            map_nonstd_res=True,
            assign_elem=True,
        )

        mol.alg.Molck(
            self.processed_model,
            conop.GetDefaultLib(),
            molck_settings,
            prune=True,
        )
        ost.PopLogSink()
        ost.PopVerbosityLevel()
        self.preprocessing_log["molck_log"] = molcklogger.full_message
        self.preprocessing_log["atomsRemoved"] = pre_atom_count - len(
            self.processed_model.atoms
        )
        self.preprocessing_log["residuesRemoved"] = pre_res_count - len(
            self.processed_model.residues
        )
        self.preprocessing_log["chainsRemoved"] = pre_ch_count - len(
            self.processed_model.chains
        )
        self.peptide_processed_model = self.processed_model.Select(
            "peptide=true"
        )

        #########################################################
        # Add chain names if not set (happens with CASP models) #
        #########################################################
        stupid_names = (
            "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_-"
        )
        ed = self.processed_model.handle.EditXCS()
        chain_name_mapping = dict()
        for ch in self.processed_model.chains:
            if len(ch.GetName().strip()) < 1:
                # iterate over possible chain names and assign first unused name
                for new_cname in stupid_names:
                    try:
                        old_cname = ch.GetName()
                        ed.RenameChain(ch, new_cname)
                        chain_name_mapping[old_cname] = new_cname
                        break
                    except:
                        continue
            else:
                # chain name is valid, nothing to do for this chain
                chain_name_mapping[ch.GetName()] = ch.GetName()
        self.preprocessing_log["chain_name_mapping"] = chain_name_mapping

        #########################
        # Etract SEQRES/ATOMSEQ #
        #########################
        alignments = list()
        if self.seqres:
            # SEQRES is provided by user, requires mapping to model chains

            # option 1: all chains align to this single SEQRES (monomer or homo-oligomer)
            if len(self.seqres) == 1:
                for ch in self.peptide_processed_model.chains:
                    try:
                        aln = AlignChainToSEQRES(ch, self.seqres[0])
                        alignments.append(aln)
                    except:
                        raise RuntimeError(
                            f"Failed to align chain {ch.GetName()} of {self.model_path} to provided SEQRES."
                        )

            # option2: Map chains using names (whatever-mer, as long as
            #          as I find a SEQRES for each chain based on name matching)
            elif len(self.seqres) > 1:
                for ch in self.peptide_processed_model.chains:
                    ch_seqres = None
                    for s in self.seqres:
                        if ch.GetName() == s.GetName():
                            ch_seqres = s
                            break
                    if ch_seqres is None:
                        raise RuntimeError(
                            f"Failed to find SEQRES for chain of name {ch.GetName()} in provided SEQRES list."
                        )
                    try:
                        aln = AlignChainToSEQRES(ch, ch_seqres)
                        alignments.append(aln)
                    except:
                        raise RuntimeError(
                            f"Failed to align chain {ch.GetName()} of {self.model_path} to provided SEQRES."
                        )
        else:
            # No SEQRES provided, extract SEQRES from protein => SEQRES==ATOMSEQ
            for ch in self.peptide_processed_model.chains:
                s = "".join([r.one_letter_code for r in ch.residues])
                aln = seq.CreateAlignment()
                aln.AddSequence(seq.CreateSequence("SEQRES", s))
                aln.AddSequence(seq.CreateSequence("ATOMSEQ", s))
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
    def __init__(
        self,
        method,
        model_paths,
        seqres,
        workdir,
        uniclust30,
        smtldir,
        datefilter,
    ):
        self.method = method
        self._seqres_list = None

        # setup models
        self.models = list()
        for mp in model_paths:
            self.models.append(ModelScorer(mp, seqres))

        # set SEQRES specific data (dicts with key: seqres name, value: data):
        # - self.psipred_handler
        # - self.accpro_handler
        # - self.disco_container
        self._do_seqres_features(workdir, uniclust30, smtldir, datefilter)

        # perform scoring on all models
        self._score()

    def to_json(self):
        out_dict = dict()
        out_dict["models"] = dict()
        for m_idx, m in enumerate(self.models):
            identifier = str(m_idx + 1)
            identifier = (3 - len(identifier)) * "0" + identifier
            out_dict["models"]["model_" + identifier] = m.to_json()
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

    def _do_seqres_features(self, workdir, uniclust30, smtldir, datefilter):
        self.psipred_handler = dict()
        self.accpro_handler = dict()
        self.disco_container = dict()
        for s in self.seqres_list:
            p, a, d = _seqanno(
                s,
                workdir,
                uniclust30,
                self.method == "QMEANDisCo",
                smtldir,
                datefilter,
            )
            self.psipred_handler[s.GetName()] = p
            self.accpro_handler[s.GetName()] = a
            self.disco_container[s.GetName()] = d

    def _score(self):
        for m in self.models:
            m.score(
                self.psipred_handler,
                self.accpro_handler,
                self.disco_container,
                self.method,
            )


def _get_uniclust30():
    # uniclust30 is expected to be mounted at /uniclust30
    # however, we don't know the prefix...
    if not os.path.exists("/uniclust30"):
        raise RuntimeError("You must mount UniClust30 to /uniclust30")

    expected_uniclust30_suffixes = [
        "_a3m.ffdata",
        "_a3m.ffindex",
        "_hhm.ffdata",
        "_hhm.ffindex",
        "_cs219.ffdata",
        "_cs219.ffindex",
    ]

    uniclust_files = os.listdir("/uniclust30")

    # go over all files, check whether they have a valid suffix
    # and basically collect for each unique prefix, all observed suffixes
    prefixes = dict()
    for f in uniclust_files:
        for s in expected_uniclust30_suffixes:
            if f.endswith(s):
                prefix = f[: -len(s)]
                if prefix not in prefixes:
                    prefixes[prefix] = list()
                prefixes[prefix].append(s)
                break

    # we can now check whether the expected suffixes are complete for
    # each unique prefix
    complete_prefixes = list()
    for prefix, suffixes in prefixes.items():
        if sorted(suffixes) == sorted(expected_uniclust30_suffixes):
            complete_prefixes.append(os.path.join("/uniclust30", prefix))

    if len(complete_prefixes) == 1:
        return complete_prefixes[0]
    elif len(complete_prefixes) == 0:
        msg = "Expect valid UniClust30 to be mounted at /uniclust30. Files "
        msg += "with equal prefix must be present for the following suffixes: "
        msg += ", ".join(expected_uniclust30_suffixes)
        raise RuntimeError(msg)
    else:
        msg = "Expect valid UniClust30 to be mounted at /uniclust30. Files "
        msg += "with equal prefix must be present for the following suffixes: "
        msg += ", ".join(expected_uniclust30_suffixes)
        msg += "\n Several prefixes fulfilled this criterium: "
        msg += ", ".join(complete_prefixes)
        raise RuntimeError(msg)


def _check_smtl(args):
    # expect smtl to be mounted at /smtl
    if not os.path.exists("/smtl"):
        raise RuntimeError(
            "For running QMEANDisCo you need to mount the downloadable SMTL data to /smtl"
        )

    expected_files = [
        "smtl_uniq_cs219.ffdata",
        "smtl_uniq_cs219.ffindex",
        "smtl_uniq_hhm.ffdata",
        "smtl_uniq_hhm.ffindex",
        "CHAINCLUSTERINDEX",
        "indexer.dat",
        "seqres_data.dat",
        "atomseq_data.dat",
        "ca_pos_data.dat",
        "VERSION",
    ]

    for f in expected_files:
        p = os.path.join("/smtl", f)
        if not os.path.exists(p):
            raise RuntimeError(f"expect {p} to be present")

    # only if date-filter is specified, we additionally need dates.csv'
    if args.datefilter:
        p = os.path.join("/smtl", "dates.csv")
        if not os.path.exists(p):
            raise RuntimeError(
                f"If datefilter argument is provided, you additionally need to provide the SMTL specific {p}"
            )


def _parse_args():

    ####################
    # Define Arguments #
    ####################

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "models",
        nargs="+",
        help="Models in PDB or mmCIF format (also compressed, i.e. mmcif.gz)",
    )
    parser.add_argument(
        "--method",
        dest="method",
        choices=["QMEAN", "QMEANDisCo", "QMEANBrane"],
        default="QMEANDisCo",
        help="Used scoring function - Default: QMEANDisCo",
    )
    parser.add_argument(
        "--out",
        dest="out",
        default="out.json",
        help="Destination for JSON formatted output",
    )
    parser.add_argument(
        "--seqres",
        dest="seqres",
        default=None,
        help="SEQRES for models in FASTA format - Single sequence for homomers/homo-oligomers - Multiple sequences for hetero-oligomers with name based matching",
    )
    parser.add_argument(
        "--profiles",
        nargs="+",
        default=None,
        help="Precomputed HHblits sequence profile(s) in a3m format that match target sequence(s) provided in seqres - must contain psipred annotation",
    )
    parser.add_argument(
        "--workdir",
        dest="workdir",
        default=None,
        help="Location for intermediate output, normally temporary. If given, output remains for debug purposes",
    )
    parser.add_argument(
        "--datefilter", dest="datefilter", default=None, help="Debug purposes"
    )
    parser.add_argument(
        "--version",
        dest="version",
        action="store_true",
        help="Display version and exit",
    )
    args = parser.parse_args()

    # if version flag is set, we just print some stuff and abort
    if args.version:
        print("QMEAN container entry point running")
        print("QMEAN:", os.getenv("VERSION_QMEAN"))
        print("OPENSTRUCTURE:", os.getenv("VERSION_OPENSTRUCTURE"))
        sys.exit(0)

    ######################################
    # Argument checks and postprocessing #
    ######################################

    if len(args.models) >= 1000:
        raise RuntimeError("Can only score a maximum of 999 models at once")

    for model_path in args.models:
        if not os.path.exists(model_path):
            raise RuntimeError(f"specified path {model_path} does not exist")

    if args.workdir:
        # user defined workdir
        # create if it doesnt exist and disable cleanup in the end
        if not os.path.exists(args.workdir):
            os.makedirs(args.workdir)
        args.cleanup_workdir = False
        ost.LogInfo(f"User defined workdir: {args.workdir}")
    else:
        # use tmp directory as workdir which is cleaned up in the end
        args.workdir = tempfile.mkdtemp()
        args.cleanup_workdir = True
        ost.LogInfo(f"Tmp workdir: {args.workdir}")

    if args.seqres:
        if not os.path.exists(args.seqres):
            raise RuntimeError(f"specified path {args.seqres} does not exist")

    # SEQRES gets already loaded, models come later
    seqres = None
    if args.seqres:
        try:
            seqres = io.LoadSequenceList(args.seqres)
        except:
            raise RuntimeError(f"failed to load seqres file: {args.seqres}")
    args.seqres = seqres

    # Deal with profiles
    if args.profiles:
        if not args.seqres:
            raise RuntimeError("Must set seqres if providing profiles")
        seqres_hashes = [_get_seq_name(str(s)) for s in args.seqres]
        for p in args.profiles:
            if not os.path.exists(p):
                raise RuntimeError(f"specified path {p} does not exist")
            a3m_content = hhblits3.ParseA3M(open(p))
            if a3m_content["ss_pred"] is None or a3m_content["ss_conf"] is None:
                raise RuntimeError(
                    f"Sequence profile {p} must contain secondary structure annotation"
                )
            trg_seq = a3m_content["msa"].GetSequence(0).GetGaplessString()
            trg_seq_hash = _get_seq_name(trg_seq)
            if trg_seq_hash not in seqres_hashes:
                raise RuntimeError(f"Could not find matching SEQRES for {p}")
            seqanno_workdir = os.path.join(args.workdir, trg_seq_hash)
            os.makedirs(seqanno_workdir)
            shutil.copy(
                p,
                os.path.join(
                    seqanno_workdir, hhblits3.HHblits.OUTPUT_PREFIX + ".a3m"
                ),
            )

    # Check for valid uniclust30 and report to logger
    args.uniclust30 = _get_uniclust30()
    ost.LogInfo(f"Use UniClust30: {args.uniclust30}")

    # Check for valid SMTL in case of QMEANDisCo
    if args.method == "QMEANDisCo":
        _check_smtl(args)

    # check what complib is used and report to logger
    creation_date = conop.GetDefaultLib().GetCreationDate()
    ost.LogInfo(f"Use compound library with creation date {creation_date}")

    return args


def _main():

    # start logging
    logger = _Logger()
    ost.PushLogSink(logger)
    ost.PushVerbosityLevel(ost.LogLevel.Info)

    ##############################################
    # Parse arguments and already prepare output #
    ##############################################
    args = _parse_args()
    out = dict()
    out["qmean_version"] = qmean.qmean_version
    out["created"] = datetime.datetime.now().isoformat(timespec="seconds")
    out["method"] = args.method
    out["seqres_uploaded"] = None
    if args.seqres:
        out["seqres_uploaded"] = list()
        for s in args.seqres:
            out["seqres_uploaded"].append(
                {"name": s.GetName(), "sequence": str(s)}
            )
    out["smtl_version"] = None
    if args.method == "QMEANDisCo":
        with open("/smtl/VERSION", "r") as fh:
            out["smtl_version"] = fh.read().strip()

    ######################################
    # Load/process models and do scoring #
    ######################################
    data = ModelScorerContainer(
        args.method,
        args.models,
        args.seqres,
        args.workdir,
        args.uniclust30,
        "/smtl",
        args.datefilter,
    )

    ###############################
    # Finalize output and dump it #
    ###############################
    out["error_log"] = logger.error_message
    out["info_log"] = logger.info_message
    data_out = data.to_json()
    out.update(data_out)
    with open(args.out, "w") as fh:
        json.dump(out, fh)

    ###########
    # Cleanup #
    ###########
    if args.cleanup_workdir:
        shutil.rmtree(args.workdir)


if __name__ == "__main__":
    _main()
