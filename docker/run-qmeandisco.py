"""usage: sm run-qmeandisco <project>
"""
import re
import argparse
import os
import sys
import json

from qmean import QMEANScorer
from qmean import AssessMembraneModelQuality
from qmean import PSIPREDHandler
from qmean import ACCPROHandler
from qmean import predicted_sequence_features
from qmean import DisCoContainer
from qmean import mqa_result_membrane

import ost
from ost import io, conop
from ost.io import LoadMMCIF
from ost.bindings import hhblits
from ost.bindings import WrappedTMAlign
from ost.seq import CreateSequenceList
from ost.mol.alg import FindMembrane, MolckSettings, Molck
from ost.table import *
from ost import LogWarning

#from sm.core import cmptcmpt
#from sm.pipeline import ProjectFactory
#from sm.pipeline import RunHHblits
#from sm.pipeline import CreateDistanceConstraints
#from sm import config
#from sm.meld.qa import FetchScoresFromScorer

def _ParseArgs():
    parser = argparse.ArgumentParser(formatter_class=\
                                     argparse.ArgumentDefaultsHelpFormatter,
                                     description=__doc__)
    parser.add_argument("model", help="Model file")
    parser.add_argument("-v", "--verbose", dest="verbose_flag",
                      default=0, help="Print more detailed output (0/1)")

    opts = parser.parse_args()

    return opts


def _GetSequenceFeatures(proj):

    seqres = dict() 
    psipred_handler = dict()
    accpro_handler = dict()
    dc_container = dict()
    
    for subproject_name, subproject in proj.subprojects.items():

        seqres[subproject_name] = subproject.target_seq

        a3m_filename = subproject.a3m_filename
        acc_filename = subproject.acc_filename
        dc_filename = subproject.qmean_dc_filename

        # try to get everything for PSIPREDHandler
        if a3m_filename:
            a3m_content = hhblits.ParseA3M(open(a3m_filename))
            data = dict()
            data['ss'] = a3m_content['ss_pred']
            data['conf'] = [int(c) for c in a3m_content['ss_conf']]
            data['seq'] = str(subproject.target_seq)
            psipred_handler[subproject_name] = PSIPREDHandler(data)
        else:
            psipred_handler[subproject_name] = None

        # try to get everything for ACCPROHandler
        if acc_filename:
            data = dict()
            data["acc"] = str(io.LoadSequence(acc_filename))
            data["seq"] = str(subproject.target_seq)
            accpro_handler[subproject_name] = ACCPROHandler(data)
        else:
            accpro_handler[subproject_name] = None

        # try to get everything for DisCoContainer
        if dc_filename:
            dc_container[subproject_name] =  \
            DisCoContainer.Load(dc_filename)
        else:
            dc_container[subproject_name] = None

    return (seqres, psipred_handler, accpro_handler, dc_container)

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


def _SuperposeOnFirstChains(models):
    # models already require chain names to be set
    
    chain_lengths = [len(m.chains[0].residues) for m in models]
    reference_idx = chain_lengths.index(max(chain_lengths))
    reference_view = models[reference_idx].CreateFullView()

    for idx in range(len(models)):
        if idx == reference_idx:
            continue
        model_view = models[idx].CreateFullView()
        try:
            tm_result = WrappedTMAlign(model_view.chains[0], reference_view.chains[0])
            models[idx].handle.EditXCS().ApplyTransform(tm_result.transform)  
        except:
            pass

    return models


def _LoadRawModels(in_dir, model_files):
    models = list()

    for m_f in model_files:
        try:
            model = LoadPDB(
                os.path.join(in_dir, m_f)).Select("peptide=true")
            models.append(model)
        except:
            models.append(None)

    return models


def _LoadProcessedModels(in_dir, model_files, seqres_lists, proj):
    models = list()

    for m_f, s_l in zip(model_files, seqres_lists):
        try:
            model = LoadPDB(os.path.join(in_dir, m_f)).Select("peptide=true")
            model = _RenumberModel(model, s_l)
            models.append(model)
        except Exception as ex:
            models.append(None)

    return_list = _SuperposeOnFirstChains(models)

    return return_list

def _AssembleResultJSON(scorer, seqres):
  """
  Given a QMEANScorer object, this function will extract global and local 
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

  :param scorer:        Scorer object from QMEAN from which to extract stuff
  :param seqres:        SEQRES sequences for every chain, this list must have
                        the exact same length as there are chains in the scored
                        model
  :type scorer:         :class:`QMEANScorer`
  :type seqres:         :class:`list` of :class:`ost.seq.SequenceHandle`

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
  except:
    global_scores["qmean6_norm_score"] = None
    global_scores["qmean6_z_score"] = None

  scores = scorer.qmean6_components
  global_scores["interaction_norm_score"] = nan_to_none(scores["interaction"])
  global_scores["interaction_z_score"] = nan_to_none(scores["interaction_z_score"])
  global_scores["cbeta_norm_score"] = nan_to_none(scores["cbeta"])
  global_scores["cbeta_z_score"] = nan_to_none(scores["cbeta_z_score"])
  global_scores["packing_norm_score"] = nan_to_none(scores["packing"])
  global_scores["packing_z_score"] = nan_to_none(scores["packing_z_score"])
  global_scores["torsion_norm_score"] = nan_to_none(scores["torsion"])
  global_scores["torsion_z_score"] = nan_to_none(scores["torsion_z_score"])
  global_scores["ss_agreement_norm_score"] = nan_to_none(scores["ss_agreement"])
  global_scores["ss_agreement_z_score"] = nan_to_none(scores["ss_agreement_z_score"])
  global_scores["acc_agreement_norm_score"] = nan_to_none(scores["acc_agreement"])
  global_scores["acc_agreement_z_score"] = nan_to_none(scores["acc_agreement_z_score"])
  global_scores["avg_local_score"] = nan_to_none(scorer.avg_local_score)
  global_scores["avg_local_score_error"] = nan_to_none(scorer.avg_local_score_error)

  local_scores = dict()
  for ch, s in zip(scorer.model.chains, seqres):
    score_list = list([None] * len(s))
    per_chain_scores = scorer.local_scores[ch.GetName()]
    for rnum, score in per_chain_scores.items():
      if rnum < 1 or rnum > len(s):
        raise RuntimeError("Observed ResNum not matching provided SEQRES!")
      score_list[rnum-1] = nan_to_none(score)
    local_scores[ch.GetName()] = score_list

  result = dict()
  result["global_scores"] = global_scores
  result["local_scores"] = local_scores
  return result


#proj = ProjectFactory.FromName(sys.argv[1])

#opts, args = _ParseArgs(sys.argv[2:])

#in_dir = proj.input_dir
#out_dir = proj.results_dir
#tmp_dir = config.QMEANSERVER_TMP_DIR

## lets read the json stuff to get all required information about the project
#project_data = proj.GetProjectJSON()

## The actual action starts here, let's put everything in
## a try block to perform appropriate error handling

#try:

#    # in case of QMEANDisCo we need a full fletched template search before
#    # obtaining the distance constraints. In case of QMEAN or QMEANBrane
#    # a sequence annotation is sufficient (aacpro and the profile from which)
#    # to read the psipred data

#    job = RunHHblits(proj, wait=True, use='PROCESS')
#    if not job.completed:
#        raise RuntimeError("HHBlits was not successful!")

#    job = CreateDistanceConstraints(proj, wait=True, use='PROCESS')
#    if not job.completed:
#        raise RuntimeError("CreateDistanceConstraints was not successful!")

#    # let's fetch the psipred handlers, accpro handlers and disco containers of
#    # the full project
#    full_seqres, full_psipred_handler, full_accpro_handler, full_dc_container =\
#    _GetSequenceFeatures(proj)

#    # lets extract the names of the models and the seqres sequences
#    model_files = list()
#    model_names = list()
#    seqres_lists = list()
#    out_paths = list()

#    for item in project_data["models"]:
#        model_files.append(str(item["modelid"]) + ".pdb")
#        model_names.append(str(item["name"]))
#        seq_list = CreateSequenceList()
#        for s in item["seqres"]:
#            seq_list.AddSequence(
#                ost.seq.CreateSequence(str(s["name"]), str(s["sequence"])))
#        seqres_lists.append(seq_list)

#    # all models get loaded prior to any processing.
#    raw_models = _LoadRawModels(in_dir, model_files)
#    processed_models = _LoadProcessedModels(in_dir, model_files, seqres_lists, proj)

#    # iterate over all models and do the quality assessment
#    for mdl_idx in range(len(model_files)):

#        model = processed_models[mdl_idx]
#        raw_model = raw_models[mdl_idx]

#        out_path = os.path.join(out_dir, model_files[mdl_idx].split('.')[0])
#        out_paths.append(out_path)
#        if not os.path.exists(out_path):
#            os.makedirs(out_path)

#        # from now on things could potentially go wrong, 
#        # so we create a per model status file
#        per_model_status_path = os.path.join(out_path, 'status')

#        # set status to running
#        _WriteModelStatus(per_model_status_path, "RUNNING")

#        try:

#            # find the model specific psipred / accpro and disco information

#            psipred_handler = list()
#            accpro_handler = list()
#            dc_container = list()

#            for seqres in seqres_lists[mdl_idx]:
#                seqres_name = None
#                for k,v in full_seqres.items():
#                    if v.GetString() == seqres.GetString():
#                        seqres_name = k
#                        break
#                if seqres_name == None:
#                    raise RuntimeError("Could not find model seqres in \
#                                        preprocessed sequences!")

#                psipred_handler.append(full_psipred_handler[seqres_name])
#                accpro_handler.append(full_accpro_handler[seqres_name])
#                dc_container.append(full_dc_container[seqres_name])

#            qmean_scorer = QMEANScorer(model, psipred=psipred_handler,
#                                       accpro=accpro_handler, dc=dc_container,
#                                       use_nn=True)
#            qmean_scorer.AssignModelBFactors()
#            score_dict = FetchScoresFromScorer(qmean_scorer, 
#                                               seqres_lists[mdl_idx],
#                                               assign_avg_local_score_error=True)
#            outfile = open(os.path.join(out_path, "score_file.json"), 'w')
#            json.dump(score_dict, outfile, indent=2, sort_keys=True)
#            outfile.close()

#            SavePDB(model, os.path.join(out_path, "model_processed.pdb"))

#            # map over bfactors from processed model to rawmodel
#            for r_processed, r_raw in zip(model.residues, raw_model.residues):
#                bfac = r_processed.atoms[0].GetBFactor()
#                for a in r_raw.atoms:
#                    a.SetBFactor(bfac)

#        except:
#            raise

#except Exception as e:                               
#    raise

def _main():
    """Run QMEAN DisCO
    """
    opts = _ParseArgs()

    mdl_ent, mdl_seq = LoadMMCIF(opts.model, seqres=True)

    mdl_ent = _RenumberModel(mdl_ent, mdl_seq)
    
    qmean_scorer = QMEANScorer(mdl_ent, psipred=None,
                               accpro=None, dc=None,
                               use_nn=True)

    result_json = _AssembleResultJSON(qmean_scorer, mdl_seq)

    print(result_json)

if __name__ == "__main__":
    _main()
