"""Run QMEAN DisCo for a mmCIF/ PDB file. Writes to the same location as where
the file comes from with '.json' as suffix.
"""
from datetime import datetime
import argparse
import json
import math
import os

from qmean import QMEANScorer
from qmean import predicted_sequence_features

import ost

# from ost import conop
from ost import (
    LogLevel,
    LogSink,
    PopLogSink,
    PopVerbosityLevel,
    PushLogSink,
    PushVerbosityLevel,
)
from ost.io import LoadMMCIF, LoadPDB, LoadSequenceList

# from ost.mol.alg import Molck, MolckSettings


def _ParseArgs():
    """Fetch arguments from the command line call.
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=__doc__,
    )
    parser.add_argument("model", help="Model file")
    parser.add_argument(
        "-f",
        "--format",
        help="Specify format of input file, mmCIF or PDB. Default is to "
        + "determine from file extension.",
        default=None,
        choices=["mmcif", "pdb"],
    )
    parser.add_argument(
        "-o", "--output", help="Path to write the JSON data to."
    )
    parser.add_argument(
        "-s",
        "--seqres",
        help="File with a list of sequences corresponding to the peptide "
        + "chains in the model.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose_flag",
        default=0,
        help="Print more detailed output (0/1)",
    )

    opts = parser.parse_args()

    return opts


def _TurnIntoQMEANServerJSON(  # pylint: disable=bad-continuation
    scores_json, ost_ent, seqres, file_nm, loaded_seqres
):
    """Add "non score" data to the QMEAN JSON output.
    """
    result_json = {
        # "2019-03-07T15:33:12.932"
        # 2020-03-06 15:53:48.925243
        "created": "{:%Y-%m-%dT%H:%M:%S.%f}".format(datetime.now()),
        "models": dict(),
        "results_page": "",
        "model_pdb": "",
        "qmean_version": os.environ.get("VERSION_QMEAN", None),
        # "model_modifications": molck_json,
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
        # scores:
        # Right now it seems to work that residue numbering from structures
        # matches SEQRES. But in the past there has been user input where that
        # was broken. So in case this happens again, here is an idea for a
        # strategy:
        # Create deep copy of entity, let the renumbered copy be equipped with
        # bfactors, run along residues, fetch bfactor from renumbered entity
        # and resnum/ chain name from deep copy to undo renumbering,
        # scorer.AssignModelBFactors() seems important for that
        "scores": scores_json,
        "chains": dict(),
    }

    for i, chn in enumerate(ost_ent.chains):
        seq = seqres[i]
        if loaded_seqres:
            aln = predicted_sequence_features.AlignChainToSEQRES(chn, str(seq))
            aseq = aln.GetSequence(1)
        else:
            aseq = seq
        mdl_for_json["chains"][seq.name] = {
            "seqres": str(seq),
            "atomseq": str(aseq),
            "name": "seq_chain_%d" % i,
        }

    if not ost_ent.GetName():
        result_json["models"]["model_001"] = mdl_for_json
    else:
        result_json["models"][ost_ent.GetName()] = mdl_for_json

    return result_json


def _RenumberModel(model, seqres_list):
    """Make residue numbers correspond to SEQRES sequence.

    We assume that the chains and the sequence list have the same order.
    """
    edi = model.handle.EditXCS()

    for chn, seqres in zip(model.chains, seqres_list):
        aln = predicted_sequence_features.AlignChainToSEQRES(chn, str(seqres))
        model_chain = model.CreateEmptyView()
        model_chain.AddChain(chn, ost.mol.INCLUDE_ALL)
        aln.AttachView(1, model_chain)
        for i in range(aln.GetLength()):
            res = aln.GetResidue(1, i)
            if res.IsValid():
                new_num = ost.mol.ResNum(i + 1)
                edi.SetResidueNumber(res.handle, new_num)

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
                   Sequences must be sorted following the order of the chains.
    :type seqres: :class:`list` of :class:`ost.seq.SequenceHandle`
    """

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
    global_scores["interaction_z_score"] = NaN_To_None(
        scores["interaction_z_score"]
    )
    global_scores["cbeta_norm_score"] = NaN_To_None(scores["cbeta"])
    global_scores["cbeta_z_score"] = NaN_To_None(scores["cbeta_z_score"])
    global_scores["packing_norm_score"] = NaN_To_None(scores["packing"])
    global_scores["packing_z_score"] = NaN_To_None(scores["packing_z_score"])
    global_scores["torsion_norm_score"] = NaN_To_None(scores["torsion"])
    global_scores["torsion_z_score"] = NaN_To_None(scores["torsion_z_score"])
    global_scores["ss_agreement_norm_score"] = NaN_To_None(
        scores["ss_agreement"]
    )
    global_scores["ss_agreement_z_score"] = NaN_To_None(
        scores["ss_agreement_z_score"]
    )
    global_scores["acc_agreement_norm_score"] = NaN_To_None(
        scores["acc_agreement"]
    )
    global_scores["acc_agreement_z_score"] = NaN_To_None(
        scores["acc_agreement_z_score"]
    )
    global_scores["avg_local_score"] = NaN_To_None(scorer.avg_local_score)
    global_scores["avg_local_score_error"] = NaN_To_None(
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
            score_list[rnum - 1] = NaN_To_None(score)
        local_scores[chn.GetName()] = score_list

    result = dict()
    result["global_scores"] = global_scores
    result["local_scores"] = local_scores
    return result


def _GetOSTEntity(modelfile, load_seqres, force):
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
        if load_seqres:
            mdl_ent, mdl_seq = LoadMMCIF(modelfile, seqres=True)
        else:
            mdl_ent = LoadMMCIF(modelfile, seqres=False)
            mdl_seq = None
    elif sformat in ["pdb"]:
        sformat = "pdb"
        if load_seqres:
            mdl_ent, mdl_seq = LoadPDB(modelfile, seqres=True)
        else:
            mdl_ent = LoadPDB(modelfile, seqres=False)
            mdl_seq = None
    else:
        raise RuntimeError(
            "Unknow/ unsuppoted file extension found for file '%s'" % modelfile
        )

    return mdl_ent.Select("peptide=true"), mdl_seq, sformat


class _MolckToJsonLogger(LogSink):
    """Append Molck output to results object."""

    def __init__(self):
        LogSink.__init__(self)
        self.processed = _MolckToJsonLogger.SetupDict()

    def __enter__(self):
        PushVerbosityLevel(LogLevel.Info)
        PushLogSink(self)
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        PopVerbosityLevel()
        PopLogSink()

    def LogMessage(self, message, severity):  # pylint: disable=unused-argument
        """Send a message to the logger.
        """
        message = message.strip()
        if message.endswith(" is not a standard amino acid --> removed"):
            res = message.split()[0]
            self.processed["removed_non_aa"].append(res)
        else:
            raise RuntimeError("Found unknown Molck output: '%s'" % message)

    def ToJson(self):
        """Return what the logger collected as JSON.
        """
        return self.processed

    @staticmethod
    def SetupDict():
        """Create the empty JSON structure.
        This is meant for the case the logger is not in use but you still want
        the keys, etc. of the JSON output for the full QMEAN JSON format.
        """
        return {"removed_non_aa": []}


def _main():
    """Run QMEAN DisCO
    """
    opts = _ParseArgs()

    load_seqres = True
    if opts.seqres:
        load_seqres = False
    mdl_ent, mdl_seq, fmt = _GetOSTEntity(opts.model, load_seqres, opts.format)

    # We need to be careful with Molck since it alters the input structure.
    # There should not be a problem with proper mmCIF files but with PDB.
    # rm_non_std - removes non-amino acid residues like HOH and ligands. If
    #              those residues appear within a peptide chain, OST gets
    #              problems dealing with the sequence corresponding to that
    #              chain. It is safe to remove water and ligands since QMEAN
    #              does not touch them.
    # molck_out = _MolckToJsonLogger.SetupDict()

    # if fmt == "pdb":
    #    with _MolckToJsonLogger() as molcklogger:
    #        Molck(
    #            mdl_ent,
    #            conop.GetDefaultLib(),
    #            MolckSettings(
    #                rm_unk_atoms=False,
    #                rm_non_std=True,
    #                rm_hyd_atoms=False,
    #                rm_oxt_atoms=False,
    #                rm_zero_occ_atoms=False,
    #                colored=False,
    #                map_nonstd_res=False,
    #                assign_elem=False,
    #            ),
    #        )
    #        molck_out = molcklogger.ToJson()

    if opts.seqres:
        mdl_seq = LoadSequenceList(opts.seqres, format="fasta")

    mdl_ent = _RenumberModel(mdl_ent, mdl_seq)

    qmean_scorer = QMEANScorer(
        mdl_ent, psipred=None, accpro=None, dc=None, use_nn=True
    )

    scores_json = _AssembleResultJSON(qmean_scorer, mdl_seq)

    result_json = _TurnIntoQMEANServerJSON(
        scores_json,
        mdl_ent,
        mdl_seq,
        os.path.basename(opts.model),
        not load_seqres,
    )

    if not opts.output:
        out_path = "%s.json" % os.path.splitext(opts.model)[0]
    else:
        out_path = opts.output

    with open(out_path, "w") as ofh:
        json.dump(result_json, ofh, indent=2, sort_keys=True)


if __name__ == "__main__":
    _main()

#  LocalWords:  OST mmCIF PDB param Molck HOH ligands QMEAN SEQRES
