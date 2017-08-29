import os, pickle, traceback
import sklearn
from sklearn.ensemble import RandomForestRegressor
from non_redundant_model import Pick
from qmean import predicted_sequence_features, mqa_result
import qmean


# This file needs to be dumped into an unpacked version of the massive 
# qmeandisco_SMset_training_and_evaluation_data.tar.gz that is stored
# in Christine Rempfers home => better data location required!
# You simply run by:
# ost tree_data_extraction.py


def DumpCSV(target_file, csv_outfile):
  """
  target_file: file with one smng project name per line (e.g. 1jcd_A)
  csv_outfile: the name of the CSV file that will be produced for training and
            testing the qmeandisco random forest
  """

  # compute features neccessary for training and evaluation of the random forest 
  # from DCData as well as QMEAN and DisCo scores
  data = list()
  targets = open(target_file,'r').readlines()
  targets = [tt.strip() for tt in targets]

  print "OBTAIN DATA FOR TARGETS IN: ", target_file
  print "NUM TARGETS: ", len(targets)

  for tt_nr, tt in enumerate(targets):
    try:
      print "PROCESSING ",tt, tt_nr

      path_to_smng_project = os.path.join("testset","smng_projects2",tt.strip()+".sm")
      target_seq = io.LoadSequence(os.path.join(path_to_smng_project,"target.fasta"))
      models = Pick(path_to_smng_project)

      disco_data = qmean.LoadDCData("dc_data/%s.dat" % tt.strip())

      print "PICKED MODELS: ",models

      infile = open(os.path.join("testset","lddt_data",tt.strip()+".dat"),'rb')
      lddt_data = pickle.load(infile)
      infile.close()

      infile = open(os.path.join("testset","ss_pred_data",tt.strip()+".dat"),'rb')
      ss_acc_data = pickle.load(infile)
      infile.close()

      psi = predicted_sequence_features.PSIPREDHandler(ss_acc_data)
      acc = predicted_sequence_features.ACCPROHandler(ss_acc_data)    

      for mdl in models:
        try:
          path_to_model = os.path.join(path_to_smng_project,"mdl",mdl,"model.pdb")
          prot = io.LoadPDB(path_to_model).Select("peptide=true and ligand=false")

          result = mqa_result.AssessModelQuality(prot,psipred=psi,accpro=acc,plots=False) 
          qmean_score = list()
          for res in prot.residues:
            qmean_score.append(res.atoms[0].b_factor)
          global_qmean4 = result[0].qmean4.norm        

          local_lddt = lddt_data[mdl]["local_lddt"]

          aln = seq.alg.AlignToSEQRES(prot, target_seq)
          aln.AttachView(1,prot)

          mdl_features = qmean.DetermineFeatureValues(aln, disco_data)
          mdl_features = [list(x) for x in mdl_features] #convert to python lists
          disco = qmean.DCScore(aln, disco_data)

          for j,res_features in enumerate(mdl_features):        
            row = [tt,mdl, res.number.GetNum(),"90",float(local_lddt[j]), global_qmean4]
            row.extend(res_features)
            row.append(qmean_score[j])
            row.append(disco[j])
            data.append(row) 

        except:
          traceback.print_exc() 

    except:
      traceback.print_exc()

  columns = ["Target","Model","ResNum","SeqIdCo","lDDT","GlobalQMEAN4","Avg.Max.SeqSim", 
             "NumClust", "Variance","NumInter", "Avg.Max.SeqId","QMEAN", "DisCo"]

  outfile = open(csv_outfile, 'w')

  outfile.write(','.join(columns))
  outfile.write('\n')

  for data_line in data:
    outfile.write(','.join([str(item) for item in data_line]))
    outfile.write('\n')

  outfile.close() 


DumpCSV("qmeandisco_training_targets.txt", "tree_training_data.csv")
DumpCSV("scorer_test_targets.txt", "tree_test_data.csv")

