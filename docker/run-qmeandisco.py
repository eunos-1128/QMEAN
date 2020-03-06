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
from ost.io import (LoadPDB,  SavePDB)
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


def _LoadMolckedModel(model_path):

    model = LoadPDB(model_path) 
    Molck(model, conop.GetDefaultLib(), MolckSettings(rm_unk_atoms=True,
                        rm_non_std=True, rm_hyd_atoms=True, rm_oxt_atoms=False,
                        rm_zero_occ_atoms=True, colored=False, map_nonstd_res=True,
                        assign_elem=True))

    return model


def _NameChainsIfNecessary(model, proj, model_file):
    valid_chain_names = True
    jquery_chars = "(?:\\\\.|[\\w-]|[^\\x00-\\xa0])+"
    for c in model.chains:
        cname = c.GetName()
        if len(cname.strip()) < 1 or re.match(jquery_chars, cname) is None:
            valid_chain_names = False
            break

    if valid_chain_names:
        return model
    ed = model.handle.EditXCS()
    name_mapping = dict()
    stupid_names = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_-"
    for c in model.handle.chains:
        for newname in stupid_names:
            try:
                name_mapping[c.GetName()] = newname
                ed.RenameChain(c, newname)
                break
            except:
                # Continue through the list until
                # we have an unused chainname
                continue
    info = proj.GetProjectJSON()
    model_id = model_file.split(".")[0]
    if "info" not in info:
        info["info"] = dict()
    if model_id not in info["info"]:
        info["info"][model_id] = dict()
    info["info"][model_id]["chain_mapping"] = name_mapping
    proj.SaveProjectJSON(info)
    name_mapping = ["%s -> %s" % (o, n) for o, n in name_mapping.items()]
    proj.UpdateProgress("renaming chains in %s: %s" % (
        model_file.split(".")[0],
        ", ".join(name_mapping)))
    return model


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
            model = _LoadMolckedModel(
                os.path.join(in_dir, m_f)).Select("peptide=true")
            models.append(model)
        except:
            models.append(None)

    return models


def _LoadProcessedModels(in_dir, model_files, seqres_lists, proj):
    models = list()

    for m_f, s_l in zip(model_files, seqres_lists):
        try:
            model = _LoadMolckedModel(
                os.path.join(in_dir, m_f)).Select("peptide=true")
            model = _NameChainsIfNecessary(model, proj, m_f)
            model = _RenumberModel(model, s_l)
            models.append(model)
        except Exception as ex:
            models.append(None)

    return_list = _SuperposeOnFirstChains(models)

    return return_list


#proj = ProjectFactory.FromName(sys.argv[1])

#opts, args = _ParseArgs(sys.argv[2:])

#in_dir = proj.input_dir
#out_dir = proj.results_dir
#tmp_dir = config.QMEANSERVER_TMP_DIR

## lets read the json stuff to get all required information about the project
#project_data = proj.GetProjectJSON()

#use_qmeandisco = project_data["options"]["qmeandisco"]

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

#        proj.UpdateProgress("processing model " + model_files[mdl_idx])
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
#                                       use_nn = use_qmeandisco)
#            qmean_scorer.AssignModelBFactors()
#            score_dict = FetchScoresFromScorer(qmean_scorer, 
#                                               seqres_lists[mdl_idx],
#                                               assign_avg_local_score_error = use_qmeandisco)
#            utfile = open(os.path.join(out_path, "score_file.json"), 'w')
#            json.dump(score_dict, outfile, indent=2, sort_keys=True)
#            outfile.close()

#            proj.UpdateProgress(" + saving down model with bfactors")
#            SavePDB(model, os.path.join(out_path, "model_processed.pdb"))

#            # map over bfactors from processed model to rawmodel
#            for r_processed, r_raw in zip(model.residues, raw_model.residues):
#                bfac = r_processed.atoms[0].GetBFactor()
#                for a in r_raw.atoms:
#                    a.SetBFactor(bfac)

#            SavePDB(raw_model, os.path.join(out_path, "model_raw.pdb"))
                
#            _WriteModelStatus(per_model_status_path, "COMPLETED")

#            #Calling LoadModelDetails will cache some alignment details
#            proj.LoadModelDetails(project_data["models"][mdl_idx]['modelid'])

#        except:
#            raise

#except Exception as e:                               
#    raise

def _main():
    """Run QMEAN DisCO
    """
    opts = _ParseArgs()

if __name__ == "__main__":
    _main()
