"""Run QMEAN DisCo for a mmCIF/ PDB file. Writes to the same location as where
the file comes from with '.json' as suffix.
"""
import argparse
import os
import json
from datetime import datetime

from qmean import QMEANScorer
from qmean import predicted_sequence_features

import ost
from ost import conop
from ost import LogLevel, LogSink
from ost.io import LoadMMCIF, LoadPDB
from ost.table import *
from ost.mol.alg import Molck, MolckSettings


def _ParseArgs():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=__doc__,
    )
    parser.add_argument("model", help="Model file")
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose_flag",
        default=0,
        help="Print more detailed output (0/1)",
    )
    parser.add_argument(
        "-o", "--output", help="Path to write the JSON data to."
    )

    opts = parser.parse_args()

    return opts


def _TurnIntoQMEANServerJSON(scores_json, ost_ent, seqres, file_nm, molck_json):
    result_json = {
        # "2019-03-07T15:33:12.932"
        # 2020-03-06 15:53:48.925243
        "created": "{:%Y-%m-%dT%H:%M:%S.%f}".format(datetime.now()),
        "models": dict(),
        "results_page": "",
        "model_pdb": "",
        "qmean_version": os.environ.get("VERSION_QMEAN", None),
        "model_modifications": molck_json,
        # Project naming may be a future feature
        "project_name": "Untitled Project",
        # This needs to be fixed once we have the data in for doing DisCo
        # sequence searches!
        "smtl_version": "2020-06-24",
        "method": "QMEAN",
        "seqres_uploaded": None,
        # Status may chainge once we have error handling
        "status": "COMPLETED",
    }

    # get per model scores
    mdl_for_json = {
        "original_name": file_nm,
        "scores": scores_json,  # scores: create deep copy of entity, let the
        # renumbered copy be equipped with bfactors, run
        # along residues, fetch bfactor from renumbered
        # entity and resnum/ chain name from deep copy
        # to undo renumbering,
        # scorer.AssignModelBFactors() seems important
        # for that
        "chains": dict(),
    }

    for i, seq in enumerate(seqres):
        mdl_for_json["chains"][seq.name] = {
            "seqres": str(seq),
            # As soon as SEQRES becomes an input of this script,
            # qmean.predicted_sequence_features.AlignChainToSEQRES needs to be
            # used. That changes the way 'atomseq' is handled. Check forms.py of
            # the QMEAN web server for this.
            "atomseq": str(seq),
            "name": "seq_chain_%d" % i,
        }

    if not ost_ent.GetName():
        result_json["models"]["model_001"] = mdl_for_json
    else:
        result_json["models"][ost_ent.GetName()] = mdl_for_json

    return result_json


def _RenumberModel(model, seqres_list):
    ed = model.handle.EditXCS()

    for ch, s in zip(model.chains, seqres_list):
        aln = predicted_sequence_features.AlignChainToSEQRES(ch, str(s))
        model_chain = model.CreateEmptyView()
        model_chain.AddChain(ch, ost.mol.INCLUDE_ALL)
        aln.AttachView(1, model_chain)
        for i in range(aln.GetLength()):
            res = aln.GetResidue(1, i)
            if res.IsValid():
                new_num = ost.mol.ResNum(i + 1)
                ed.SetResidueNumber(res.handle, new_num)

    return model


def _AssembleResultJSON(scorer, seqres):
    """ Given a QMEANScorer object, this function will extract global and local
    scores and return them as a dictionary. You can access them by the keys
    "global_scores" and "local_scores".

    For "global_scores" you get another dictionary with keys:

    * "qmean4_norm_score"
    * "qmean4_z_score"
    * "qmean6_norm_score"
    * "qmean6_z_score"
    * "interaction_norm_score"
    * "interaction_z_score"
    * "cbeta_norm_score"
    * "cbeta_z_score"
    * "packing_norm_score"
    * "packing_z_score"
    * "torsion_norm_score"
    * "torsion_z_score"
    * "ss_agreement_norm_score"
    * "ss_agreement_z_score"
    * "acc_agreement_norm_score"
    * "acc_agreement_z_score"
    * "avg_local_score"

    For "local_scores" you get another dictionary with chain names as keys
    and lists with local scores as value. The local score lists have the
    same length as the according SEQRES. Assuming that the residue numbers
    in the model correspond to their location in the SEQRES, they will determine
    the location of a per residue score in this result list. The underlying
    assumption there is a SEQRES numbering scheme starting from one.

    :param scorer: Scorer object from QMEAN from which to extract stuff.
    :type scorer: :class:`QMEANScorer`
    :param seqres: SEQRES sequences for every chain, this list must have the
                   exact same length as there are chains in the scored model.
    :type seqres: :class:`list` of :class:`ost.seq.SequenceHandle`
    """

    # QMEAN gives NaN if something is not defined. JSON prefers None
    def nan_to_none(val):
        if val != val:
            return None
        return val

    global_scores = dict()
    global_scores["qmean4_norm_score"] = nan_to_none(scorer.qmean4_score)
    global_scores["qmean4_z_score"] = nan_to_none(scorer.qmean4_z_score)

    try:
        global_scores["qmean6_norm_score"] = nan_to_none(scorer.qmean6_score)
        global_scores["qmean6_z_score"] = nan_to_none(scorer.qmean6_z_score)
    except Exception:
        global_scores["qmean6_norm_score"] = None
        global_scores["qmean6_z_score"] = None

    scores = scorer.qmean6_components
    global_scores["interaction_norm_score"] = nan_to_none(scores["interaction"])
    global_scores["interaction_z_score"] = nan_to_none(
        scores["interaction_z_score"]
    )
    global_scores["cbeta_norm_score"] = nan_to_none(scores["cbeta"])
    global_scores["cbeta_z_score"] = nan_to_none(scores["cbeta_z_score"])
    global_scores["packing_norm_score"] = nan_to_none(scores["packing"])
    global_scores["packing_z_score"] = nan_to_none(scores["packing_z_score"])
    global_scores["torsion_norm_score"] = nan_to_none(scores["torsion"])
    global_scores["torsion_z_score"] = nan_to_none(scores["torsion_z_score"])
    global_scores["ss_agreement_norm_score"] = nan_to_none(
        scores["ss_agreement"]
    )
    global_scores["ss_agreement_z_score"] = nan_to_none(
        scores["ss_agreement_z_score"]
    )
    global_scores["acc_agreement_norm_score"] = nan_to_none(
        scores["acc_agreement"]
    )
    global_scores["acc_agreement_z_score"] = nan_to_none(
        scores["acc_agreement_z_score"]
    )
    global_scores["avg_local_score"] = nan_to_none(scorer.avg_local_score)
    global_scores["avg_local_score_error"] = nan_to_none(
        scorer.avg_local_score_error
    )

    local_scores = dict()
    for chn, seq in zip(scorer.model.chains, seqres):
        score_list = list([None] * len(seq))
        per_chain_scores = scorer.local_scores[chn.GetName()]
        for rnum, score in per_chain_scores.items():
            if rnum < 1 or rnum > len(seq):
                raise RuntimeError(
                    "Observed ResNum not matching provided SEQRES!"
                )
            score_list[rnum - 1] = nan_to_none(score)
        local_scores[chn.GetName()] = score_list

    result = dict()
    result["global_scores"] = global_scores
    result["local_scores"] = local_scores
    return result


def _GetOSTEntity(modelfile, force=None):
    """Read an OST entity either from mmCIF or PDB. Also fetch the sequence.

    :param modelfile: The path to the file containing the structure.
    :type modelfile: :class:`str`
    :param force: Enforce a file format, takes "mmcif" or "pdb".
    :type force: :class:`str`
    """
    # Determine file format. If not forced, look at suffix.
    sformat = force
    if not sformat:
        ext = modelfile.split(".")
        if ext[-1] == "gz":
            ext = ext[:-1]
        if len(ext) > 1:
            sformat = ext[-1].lower()
        else:
            raise RuntimeError(
                "Could not determine format of file '%s'." % modelfile
            )

    # Load the structure
    if sformat in ["mmcif", "cif"]:
        sformat = "mmcif"
        mdl_ent, mdl_seq = LoadMMCIF(modelfile, seqres=True)
    elif sformat in ["pdb"]:
        sformat = "pdb"
        mdl_ent, mdl_seq = LoadPDB(modelfile, seqres=True)
    else:
        raise RuntimeError(
            "Unknow/ unsuppoted file extension found for file '%s'" % modelfile
        )

    return mdl_ent, mdl_seq, sformat


class _MolckToJsonLogger(LogSink):
    """Append Molck output to results object."""

    def __init__(self):
        LogSink.__init__(self)
        self.processed = _MolckToJsonLogger._SetupDict()

    def __enter__(self):
        PushVerbosityLevel(LogLevel.Info)
        PushLogSink(self)
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        PopVerbosityLevel()
        PopLogSink()

    def LogMessage(self, message, severity):
        message = message.strip()
        if message.endswith(" is not a standard amino acid --> removed"):
            res = message.split()[0]
            self.processed["removed_non_aa"].append(res)
        else:
            raise RuntimeError("Found unknown Molck output: '%s'" % message)

    def ToJson(self):
        return self.processed

    @staticmethod
    def _SetupDict():
        return {"removed_non_aa": []}


def _main():
    """Run QMEAN DisCO
    """
    opts = _ParseArgs()

    mdl_ent, mdl_seq, fmt = _GetOSTEntity(opts.model)

    # We need to be careful with Molck since it alters the input structure.
    # There should not be a problem with proper mmCIF files but with PDB.
    # rm_non_std - removes non-amino acid residues like HOH and ligands. If
    #              those residues appear within a peptide chain, OST gets
    #              problems dealing with the sequence corresponding to that
    #              chain. It is safe to remove water and ligands since QMEAN
    #              does not touch them.
    molck_out = _MolckToJsonLogger._SetupDict()
    if fmt == "pdb":
        with _MolckToJsonLogger() as molcklogger:
            Molck(
                mdl_ent,
                conop.GetDefaultLib(),
                MolckSettings(
                    rm_unk_atoms=False,
                    rm_non_std=True,
                    rm_hyd_atoms=False,
                    rm_oxt_atoms=False,
                    rm_zero_occ_atoms=False,
                    colored=False,
                    map_nonstd_res=False,
                    assign_elem=False,
                ),
            )
            molck_out = molcklogger.ToJson()

    mdl_ent = _RenumberModel(mdl_ent, mdl_seq)

    qmean_scorer = QMEANScorer(
        mdl_ent, psipred=None, accpro=None, dc=None, use_nn=True
    )

    scores_json = _AssembleResultJSON(qmean_scorer, mdl_seq)

    result_json = _TurnIntoQMEANServerJSON(
        scores_json, mdl_ent, mdl_seq, os.path.basename(opts.model), molck_out
    )

    if not opts.output:
        out_path = "%s.json" % os.path.splitext(opts.model)[0]
    else:
        out_path = opts.output

    with open(out_path, "w") as ofh:
        json.dump(result_json, ofh, indent=2, sort_keys=True)


if __name__ == "__main__":
    _main()

#  LocalWords:  OST mmCIF PDB param Molck HOH ligands QMEAN
