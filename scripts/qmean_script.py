#!/usr/bin/env ost
import os, os.path, sys
from optparse import OptionParser,SUPPRESS_HELP

from ost import io,settings

from ost.io import repository
import qmean
from qmean import statistical_potentials
from qmean import conf, global_scores, local_scores
from qmean import predicted_sequence_features
#import cProfile

"""
Implementation of QMEAN in OST. Provide one or several PDBs and the QMEAN score is calculated
and compared to X-ray structure of similar size. A plot is generated in R visualising the relationship.
A Z-score of the QMEAN QMEAN-score (and all contributing terms) of the structure with respect to the
reference X-rays is provided.

Usage:
ost qmean_script.py -h

Author: Gabriel Studer, gabriel.studer@stud.unibas.ch
"""

def _ParseOptions():
  parser=OptionParser()
  parser.add_option('-l', '--pdb_list', default=None, type='string', 
                    help='list of PDB ids (e.g. PISCES set or set of PDB+chain)')
  parser.add_option('-d', '--pdb_dir', default=None, type='string', 
                    help='directory containing PDB-files')
  parser.add_option('-p', '--pdb_path', default=None, help='Path to a PDB of a PDB identifier')                
  parser.add_option('-c', '--chain', type='string', default='')
  parser.add_option('-o', '--output_dir', default=".")
  parser.add_option('-I', '--identifier', default="Txxxx", help='identifier used as basename for all output files')
  parser.add_option('--dssp-bin', default=None, type='string', help='path to dssp binary')
  #agreement terms
  parser.add_option('-S', '--psipred', default=None, help='PSIPRED file (currently only fasta format supported)')
  parser.add_option('-A', '--sspro', default=None, help='SSpro/ACCpro file (*.sspro_acc)')
  parser.add_option('-a', '--predicted_features_directory', default=None, help='Directory with PSIPRED and SSpro/ACCpro files (MD5.SSpro_acc, MD5.horiz)') 

  parser.add_option('-f', action="store_true", dest="use_only_first_chain", default=False, help='use only first chain, ignore others')
  parser.add_option('-C', action="store_true", dest="analyse_entire_complex", default=False, help='analyse protein complex as a whole, including all chains')
  parser.add_option('--additional_chain_scores', action="store_true", default=False, help='for protein complex, also scores for individual chains are delivered')

  parser.add_option('-L', action="store_true", dest='calculate_local_scores', help='Calculate local scores')
  parser.add_option("-w", "--window_size", default=9, help=SUPPRESS_HELP)
  parser.add_option("--spherical_smoothing_mode", action="store_true", default=False, help=SUPPRESS_HELP)
  parser.add_option("-W", "--normalisation_flag", action="store_true", default=False, help=SUPPRESS_HELP)
  parser.add_option("-R", '--quality_R_plots', default=0, help='Quality of Z-score plots (generated in R): 0 png,  1 Cairo (w/o X11), 2 tiff')
  parser.add_option('-P', action="store_true", dest='plot_energy_profile', help='Plot predicted local model quality ("energy profile")')
  parser.add_option('-B', action="store_true", dest='save_Bfactor_colored_PDB', help='Color the model according to local error ("energy profile") in B-factor column')
  parser.add_option("-Z", action="store_false", dest='slider_flag', help="Skipp generation of Z-score slider file", default=True)
  parser.add_option("-D", '--local_distance_measure', default=1, help='Distance measure for local error: 0 S-score,  1 RMSD (default)')
  parser.add_option('-X', action="store_true", dest='CASP_output', help=SUPPRESS_HELP)
  parser.add_option('-Q', '--method', default="QMEAN", help=SUPPRESS_HELP)

  parser.add_option("-m", action="store_true", dest="merge_output", default=False, help='show multiple structures in one plot')
  parser.add_option("-z", action="store_true", dest="label_structures", default=False, help='label structures in Z-score plots')
  parser.add_option("-n", action="store_true", dest="no_plots", default=False, help='no Z-score plots are generated')
  parser.add_option("-t", action="store_true", dest="save_thumbnail", default=False, help='create thumbnail 150x150')

  parser.add_option("-V", action="store_true", dest="verbose_output_mode", default=False, help='include full data in output table')
  parser.add_option("-N", action="store_false", dest="list_has_header", default=True, help='for -l option: list file contains no header')
  parser.add_option("-s", action="store_true", dest="silent_mode", default=False, help='do not show any scores on stdout')

  parser.add_option("-M", action="store_true", dest="Membrane_potentials", default=False, help=SUPPRESS_HELP)
                    #help='use potentials trained on membrane proteins')
  parser.add_option("-x", action="store_true", dest="assign_experimental_method", default=False, help=SUPPRESS_HELP)
                    #help='assign experimental method derived from PDB/derived_data')
  parser.add_option("-F", action="store_true", dest="fraction_modelled_flag", default=False, help=SUPPRESS_HELP)
                    #help='Divide QMEAN scores by fraction modelled (QMEANserver mode only)')
  parser.add_option("-O", action="store_true", dest="use_occupancy_flag", default=False, help=SUPPRESS_HELP)
                    #help='Take into account occupancy column in PDB file')
  parser.add_option("-r", "--repository", default="PDB", help=SUPPRESS_HELP)
                    #help="interaction mode: PDB, EXPDB, PDBredo or PISA [default: %default]")
  parser.add_option("--no_zscores", default=False, action="store_true", help=SUPPRESS_HELP)
  parser.add_option("--QMEAN5", action="store_true", default=False, help=SUPPRESS_HELP) #no all-atom term



  opts, args=parser.parse_args()

  if len(sys.argv) < 2:
    print "\nSpecify file(s) -> option -h for help\n"
    sys.exit(0)

  if opts.output_dir != ".":
    os.system("mkdir -p "+opts.output_dir)

  return opts




def _initialiseRepository(repo):

  if 'BC2_DATA_DIR' not in globals():
    raise RuntimeError("No access to PDB repository. On login define 'BC2_DATA_DIR' in .ostrc")

  import string
  #if opts.silent_mode == False: 
    #print "repository:",opts.repository

  #repository for original (all-chain) PDBs:
  if opts.analyse_entire_complex== True:
    if repo=="PISA":
      print "PISA not supported yet"
    elif repo == "EXPDB":
      #print "For using the EXPDB, a chain identifier needs to be defined. PDB taken instead."
      PDB_repo=repository.ModelRepository(os.path.join(BC2_DATA_DIR,'PDB/data/structures/all/pdb'), 
                                  file_pattern='pdb%(id)s.ent.gz', 
                                  transform=string.lower)
    elif repo=="PDBredo":
      PDB_repo=repository.ModelRepository(os.path.join(BC2_DATA_DIR,'PDB_REDO/'), 
                               file_pattern='%(dir)s/%(id)s/%(id)s_besttls.pdb', 
                               transform=string.lower)

    #elif repo=="biounit":
      #PDB_repo=repository.ModelRepository(os.path.join(BC2_DATA_DIR,'PDB/data/biounit/coordinates/all'), 
                               #file_pattern='%(id)s.ent.gz', 
                               #transform=string.lower)

    else:
      PDB_repo=repository.ModelRepository(os.path.join(BC2_DATA_DIR,'PDB/data/structures/all/pdb'), 
                                  file_pattern='pdb%(id)s.ent.gz', 
                                  transform=string.lower)
  else:
    if repo=="PDB":
      PDB_repo=repository.ModelRepository(os.path.join(BC2_DATA_DIR,'PDB/data/structures/all/pdb'), 
                                  file_pattern='pdb%(id)s.ent.gz', 
                                  transform=string.lower)
    elif repo=="PISA":
      print "PISA not supported yet"

    elif repo=="PDBredo":
      PDB_repo=repository.ModelRepository(os.path.join(BC2_DATA_DIR,'PDB_REDO/'), 
                               file_pattern='%(dir)s/%(id)s/%(id)s_besttls.pdb', 
                               transform=string.lower)
    else: #EXPDB
      #if opts.chain == '' and opts.pdb_list == None and opts.pdb_dir == None:        #for list files chains are not defined
      if opts.chain == '':        #for list files chains are not defined
        #print "For using the EXPDB, a chain identifier needs to be defined. PDB taken instead."
        PDB_repo=repository.ModelRepository(os.path.join(BC2_DATA_DIR,'PDB/data/structures/all/pdb'), 
                                  file_pattern='pdb%(id)s.ent.gz', 
                                  transform=string.lower)
      else:
        PDB_repo=repository.ModelRepository(os.path.join(BC2_DATA_DIR,'EXPDB/EXPDB/'), 
                                  file_pattern='%(id)s%(chain)s.pdb', 
                                  transform=string.lower)
  return PDB_repo



def _SinglePdb(pdb_id, chain):
  global qmean_scores
  clean_up_flag=0
  if len(pdb_id) == 4:
    path=pdb_id
    PDB_repo=_initialiseRepository(opts.repository)
    ent=PDB_repo.Load(pdb_id, chain)
  else:
    path=pdb_id
    pdb_id=os.path.basename(pdb_id)#.split('.')[0]
    ent=io.LoadPDB(str(path))
  if ent.IsValid() == 0:
    raise IOError("PDB file %s does not exist" %(pdb_id))

  out_file=os.path.join(opts.output_dir,"QMEAN_scores_"+pdb_id+".csv")
  out_global=open(out_file, 'w')

  if opts.verbose_output_mode == False:
    out_global.write('PDB chain size fraction_modelled C_beta_pairwise all_atom_pairwise solvation combined_torsion C_beta_pairwise_norm all_atom_pairwise_norm solvation_norm combined_torsion_norm SSE_agree ACC_agree QMEAN_score QMEANnorm_score\n')
  else:
    out_global.write('PDB chain size fraction_modelled fraction_residues_with_sidechains C_beta_pairwise all_atom_pairwise solvation combined_torsion C_beta_pairwise_norm all_atom_pairwise_norm solvation_norm combined_torsion_norm SSE_agree ACC_agree QMEAN4 QMEANnorm QMEAN6_score QMEANnorm6_score QMEANnorm6_score_trained_on_SwissModel QMEANnorm_size_dependent QMEANnorm6_size_dependent\n')

  #add Cbetas if missing:
  qmean.ConstructCBetas(ent, True)

  if opts.use_occupancy_flag == True:
    ent=ent.Select("occ!=0")
  SEL='peptide=true'
  ent=ent.Select(SEL)

  #if the first chain should be taken, it must be the one fitting to the SSE/ACC predition files
  #check by comparing sequence length
  if len(ent.chains) > 1 and opts.psipred!=None and opts.sspro!=None:
    y,x,pred_sequence=predicted_sequence_features.ParsePSIPRED(opts.psipred)
    pred_sequence_length=len(pred_sequence)
    for c in ent.chains:
      if c.GetResidueCount()==pred_sequence_length:
        chain=c.name
        print "Chain matching to predicted feature files:",chain
        break


  #1) scores for individual chains in multi-chain protein
  if len(ent.chains) > 1 and chain=='' and opts.analyse_entire_complex==False:
    if opts.silent_mode == False:
      #print "%s chains found..." %(len(ent.chains))
      if opts.verbose_output_mode == False:
        sys.stdout.write('PDB chain size fraction_modelled C_beta_pairwise all_atom_pairwise solvation combined_torsion C_beta_pairwise_norm all_atom_pairwise_norm solvation_norm combined_torsion_norm SSE_agree ACC_agree QMEAN_score QMEANnorm_score\n')
      else:
        sys.stdout.write('PDB chain size fraction_modelled fraction_residues_with_sidechains C_beta_pairwise all_atom_pairwise solvation combined_torsion C_beta_pairwise_norm all_atom_pairwise_norm solvation_norm combined_torsion_norm SSE_agree ACC_agree QMEAN4 QMEANnorm QMEAN6_score QMEANnorm6_score QMEANnorm6_score_trained_on_SwissModel QMEANnorm_size_dependent QMEANnorm6_size_dependent\n')
    for c in ent.chains:

      if c.GetResidueCount() == 0:
        print "Warning: PDB contains no amino acids (DNA?)"
        continue


      if opts.additional_chain_scores==True:    #analyse chain in context of entire complex
        if (opts.psipred!=None and opts.sspro!=None):
          qmean_scores=global_scores.GlobalScores((c,ent), opts.psipred, opts.sspro, stat_pot, opts.QMEAN5, dssp_path=opts.dssp_bin)
        elif opts.predicted_features_directory!=None:
          qmean_scores=global_scores.GlobalScores((c,ent), opts.predicted_features_directory, opts.predicted_features_directory, stat_pot, opts.QMEAN5, dssp_path=opts.dssp_bin)
        else:
          qmean_scores=global_scores.GlobalScores((c,ent),None,None,stat_pot, dssp_path=opts.dssp_bin)
      else: 
        if (opts.psipred!=None and opts.sspro!=None):
          qmean_scores=global_scores.GlobalScores(c, opts.psipred, opts.sspro, stat_pot, opts.QMEAN5, dssp_path=opts.dssp_bin)
        elif opts.predicted_features_directory!=None:
          qmean_scores=global_scores.GlobalScores(c, opts.predicted_features_directory, opts.predicted_features_directory, stat_pot, opts.QMEAN5, dssp_path=opts.dssp)
        else:
          qmean_scores=global_scores.GlobalScores(c,None,None,stat_pot, dssp_path=opts.dssp)

      #calculate Z-score and generate plot:
      if opts.merge_output == False:
        if  opts.use_only_first_chain==True:
          Zscore_path=os.path.join(opts.output_dir,pdb_id+"_Zscore.csv")
          plot_path=os.path.join(opts.output_dir,pdb_id+"_plot.png")
        else:
          Zscore_path=os.path.join(opts.output_dir,pdb_id+c.name+"_Zscore.csv")
          plot_path=os.path.join(opts.output_dir,pdb_id+c.name+"_plot.png")
        _calculate_Zscore_single(Zscore_path, plot_path, pdb_id+c.name)

      if opts.calculate_local_scores==True:
            predicted_features_flag=1
            if opts.psipred==None and opts.sspro==None and opts.predicted_features_directory==None:
              print "Warning: Local score calculation unreliable without PSIPRED and SSpro information"
              predicted_features_flag=0
            try:
              QMEANlocal=local_scores.LocalScores(ent, predicted_features_flag, opts.identifier, pdb_id, opts.window_size, opts.normalisation_flag ,opts.spherical_smoothing_mode, False, opts.QMEAN5)
              if opts.CASP_output==True:
                local_score_string=QMEANlocal.GetCaspVector(qmean_scores.sequence_length,1)
              #local scores output
              path=""
              model_name=""
              if  opts.use_only_first_chain==True:
                model_name=pdb_id
              else:
                model_name=pdb_id+"_"+c.name
                if c.name.strip()=="":
                  model_name_c=pdb_id+"__"
              path=os.path.join(opts.output_dir,"energy_profiles/"+model_name+"_local_energy_profile_all_terms.csv")
              if os.path.exists(path):
                os.remove(path)
              local_out_file=os.open(path, os.O_CREAT|os.O_RDWR)
              QMEANlocal.WriteOutput(local_out_file, opts.local_distance_measure)  #0:S-score, 1:RMSD
              path_PDB=os.path.join(opts.output_dir,"colored_PDB/"+model_name+"_QMEANlocal_in_Bfactor_column.pdb")
              if opts.save_Bfactor_colored_PDB==True:
                QMEANlocal.ColorPDB(path_PDB)
              os.close(local_out_file)
              if opts.plot_energy_profile==True:
                _energy_profile(path,os.path.join(opts.output_dir,"energy_profile_plots/"+model_name+"_local_energy_profile_QMEANlocal.png"), model_name)
            except Exception, e:
              print "Error in local score calculation:", e

      _chain=c.name
      if _chain=='' or _chain==" ":
          _chain="_"

      fraction_modelled=1.0
      if opts.fraction_modelled_flag == True:
        fraction_modelled=qmean_scores.fraction_modelled
      if opts.verbose_output_mode == False:
        out_global.write('%s %s %i %.2f %f %f %f %f %f %f %f %f %s %s %s %s\n' % (pdb_id,_chain,qmean_scores.protein_size,fraction_modelled, qmean_scores.energy_Cbeta, qmean_scores.energy_all_atom, qmean_scores.energy_solvation, qmean_scores.energy_torsion,qmean_scores.energy_Cbeta_norm,qmean_scores.energy_all_atom_norm,qmean_scores.energy_solvation_norm,qmean_scores.energy_torsion_norm , qmean_scores.SSE_agree, qmean_scores.ACC_agree, qmean_scores.QMEAN6_score*fraction_modelled,qmean_scores.QMEANnorm6_score*fraction_modelled))
        if opts.silent_mode == False:
          sys.stdout.write('%s %s %i %.2f %f %f %f %f %f %f %f %f %s %s %s %s\n' % (pdb_id,_chain,qmean_scores.protein_size,fraction_modelled, qmean_scores.energy_Cbeta, qmean_scores.energy_all_atom, qmean_scores.energy_solvation, qmean_scores.energy_torsion,qmean_scores.energy_Cbeta_norm,qmean_scores.energy_all_atom_norm,qmean_scores.energy_solvation_norm,qmean_scores.energy_torsion_norm , qmean_scores.SSE_agree, qmean_scores.ACC_agree, qmean_scores.QMEAN6_score*fraction_modelled,qmean_scores.QMEANnorm6_score*fraction_modelled))
      else:
        out_global.write('%s %s %i %f %f %f %f %f %f %f %f %f %f %s %s %f %f %s %s %s %f %s\n' % (pdb_id,_chain, qmean_scores.protein_size,qmean_scores.fraction_modelled,qmean_scores.fraction_residues_with_sidechains, qmean_scores.energy_Cbeta, qmean_scores.energy_all_atom, qmean_scores.energy_solvation, qmean_scores.energy_torsion,qmean_scores.energy_Cbeta_norm,qmean_scores.energy_all_atom_norm,qmean_scores.energy_solvation_norm,qmean_scores.energy_torsion_norm , qmean_scores.SSE_agree, qmean_scores.ACC_agree, qmean_scores.QMEAN_score*fraction_modelled, qmean_scores.QMEANnorm_score*fraction_modelled,qmean_scores.QMEAN6_score*fraction_modelled,qmean_scores.QMEANnorm6_score*fraction_modelled,qmean_scores.QMEANnorm6_score_trained_on_SwissModel, qmean_scores.QMEANnorm_size_dependent_score,qmean_scores.QMEANnorm6_size_dependent_score))
        if opts.silent_mode == False:
          sys.stdout.write('%s %s %i %f %f %f %f %f %f %f %f %f %f %s %s %f %f %s %s %s %f %s\n' % (pdb_id,_chain, qmean_scores.protein_size,qmean_scores.fraction_modelled,qmean_scores.fraction_residues_with_sidechains, qmean_scores.energy_Cbeta, qmean_scores.energy_all_atom, qmean_scores.energy_solvation, qmean_scores.energy_torsion,qmean_scores.energy_Cbeta_norm,qmean_scores.energy_all_atom_norm,qmean_scores.energy_solvation_norm,qmean_scores.energy_torsion_norm , qmean_scores.SSE_agree, qmean_scores.ACC_agree, qmean_scores.QMEAN_score*fraction_modelled, qmean_scores.QMEANnorm_score*fraction_modelled,qmean_scores.QMEAN6_score*fraction_modelled,qmean_scores.QMEANnorm6_score*fraction_modelled,qmean_scores.QMEANnorm6_score_trained_on_SwissModel, qmean_scores.QMEANnorm_size_dependent_score,qmean_scores.QMEANnorm6_size_dependent_score))
        
      if  opts.use_only_first_chain==True:
        break
    out_global.close()

  #1) scores for single chain or complexes
  else:
    if opts.analyse_entire_complex==False:
      full_ent=ent
      #check if chain has been specified:
      _chain=""
      if chain != '':
        _chain=chain
        ent=ent.Select("cname='%s'" %(_chain))

      if (ent.IsValid() == 0):
        print "error: chain %s not present in %s. Skipped..." %(chain, pdb_id)
        sys.exit(0)

    if opts.analyse_entire_complex== True:
      _chain=''
      # complex option renders merging useless:
      opts.merge_output=False
    if opts.additional_chain_scores==True:
      if (opts.psipred!=None and opts.sspro!=None):
        qmean_scores=global_scores.GlobalScores((ent,full_ent), opts.psipred, opts.sspro, stat_pot, opts.QMEAN5, dssp_path=opts.dssp_bin)
      elif opts.predicted_features_directory!=None:
        qmean_scores=global_scores.GlobalScores((ent,full_ent), opts.predicted_features_directory, opts.predicted_features_directory, stat_pot, opts.QMEAN5, dssp_path=opts.dssp_bin)   
      else:
        qmean_scores=global_scores.GlobalScores((ent,full_ent),None,None,stat_pot, dssp_path=opts.dssp_bin)
    else:
      if (opts.psipred!=None and opts.sspro!=None):
        qmean_scores=global_scores.GlobalScores(ent, opts.psipred, opts.sspro, stat_pot, opts.QMEAN5, dssp_path=opts.dssp_bin)
      elif opts.predicted_features_directory!=None:
        qmean_scores=global_scores.GlobalScores(ent, opts.predicted_features_directory, opts.predicted_features_directory, stat_pot, opts.QMEAN5, dssp_path=opts.dssp_bin)   
      else:
        qmean_scores=global_scores.GlobalScores(ent,None,None,stat_pot, dssp_path=opts.dssp_bin)
      #raise RuntimeError("For multichain proteins, please specify a directory containing PSIPRED and SSpro prictions in the format MD5.horiz and MD5.sspro_acc")

    #calculate Z-score and generate plot:
    Zscore_path=os.path.join(opts.output_dir,pdb_id+"_Zscore.csv")
    plot_path=os.path.join(opts.output_dir,pdb_id+"_plot.png")
    _calculate_Zscore_single(Zscore_path, plot_path, pdb_id)
    if opts.calculate_local_scores==True:
      predicted_features_flag=1
      if opts.psipred==None and opts.sspro==None and opts.predicted_features_directory==None:
        print "Warning: Local score calculation unreliable without PSIPRED and SSpro information"
        predicted_features_flag=0
      try:
        QMEANlocal=local_scores.LocalScores(ent, predicted_features_flag, opts.identifier, pdb_id, opts.window_size, opts.normalisation_flag, opts.spherical_smoothing_mode, True, opts.QMEAN5)
        if opts.CASP_output==True:
          local_score_string=QMEANlocal.GetCaspVector(qmean_scores.sequence_length)

        #local scores output
        path=""
        model_name=""
        if  opts.use_only_first_chain==True or opts.analyse_entire_complex==True:
          model_name=pdb_id
        else:
          model_name=pdb_id+_chain
        #chain_f=open(os.path.join(opts.output_dir,"chain_list.txt"),"w")
        if opts.analyse_entire_complex== True:
          chain_list=[]
          file_list=[]
          for c in ent.chains:
            model_name_c=pdb_id+"_"+c.name
            if c.name==" " or len(ent.chains)==1:
              model_name_c=pdb_id+"__"
            chain_list.append(c.name)
            #chain_f.write(c.name+";")
            path=os.path.join(opts.output_dir,"energy_profiles/"+model_name_c+"_local_energy_profile_all_terms.csv")
            file_list.append(path)
          path_PDB=os.path.join(opts.output_dir,"colored_PDB/"+model_name+"_QMEANlocal_in_Bfactor_column.pdb")
          QMEANlocal.WriteOutputComplexes(chain_list,file_list,1)  #0:S-score, 1:RMSD
          if opts.save_Bfactor_colored_PDB==True:
            QMEANlocal.ColorPDB(path_PDB)

          for i in range(len(chain_list)):
            model_name=pdb_id+"_"+chain_list[i]
            if chain_list[i]==" " or len(chain_list)==1:
              model_name=pdb_id+"__"
            if opts.plot_energy_profile==True:
              _energy_profile(file_list[i],os.path.join(opts.output_dir,"energy_profile_plots/"+model_name+"_local_energy_profile_QMEANlocal.png"), model_name)


        else:
          path=os.path.join(opts.output_dir,"energy_profiles/"+model_name+"_local_energy_profile_all_terms.csv")
          path_PDB=os.path.join(opts.output_dir,"colored_PDB/"+model_name+"_QMEANlocal_in_Bfactor_column.pdb")
          if os.path.exists(path):
            os.remove(path)
          local_out_file=os.open(path, os.O_CREAT|os.O_RDWR)
          QMEANlocal.WriteOutput(local_out_file,1)  #0:S-score, 1:RMSD
          if opts.save_Bfactor_colored_PDB==True:
            QMEANlocal.ColorPDB(path_PDB)
          os.close(local_out_file)
          if opts.plot_energy_profile==True:
            _energy_profile(path,os.path.join(opts.output_dir,"energy_profile_plots/"+model_name+"_local_energy_profile_QMEANlocal.png"), model_name)
      except Exception, e:
        print "Error in local score calculation:", e

    # for complexes -> set chain name to # chains
    if opts.analyse_entire_complex==True:
      _chain=str(ent.GetChainCount())
    if _chain=='':
      _chain="_"

    fraction_modelled=1.0
    if opts.fraction_modelled_flag == True:
      fraction_modelled=qmean_scores.fraction_modelled

    if opts.verbose_output_mode == False:
      out_global.write('%s %s %i %.2f %f %f %f %f %f %f %f %f %s %s %s %s\n' % (pdb_id,_chain, qmean_scores.protein_size,fraction_modelled, qmean_scores.energy_Cbeta, qmean_scores.energy_all_atom, qmean_scores.energy_solvation, qmean_scores.energy_torsion,qmean_scores.energy_Cbeta_norm,qmean_scores.energy_all_atom_norm,qmean_scores.energy_solvation_norm,qmean_scores.energy_torsion_norm , qmean_scores.SSE_agree, qmean_scores.ACC_agree, qmean_scores.QMEAN6_score*fraction_modelled,qmean_scores.QMEANnorm6_score*fraction_modelled))
      if opts.silent_mode == False:
        sys.stdout.write('PDB chain size fraction_modelled C_beta_pairwise all_atom_pairwise solvation combined_torsion C_beta_pairwise_norm all_atom_pairwise_norm solvation_norm combined_torsion_norm SSE_agree ACC_agree QMEAN_score QMEANnorm_score\n')
        sys.stdout.write('%s %s %i %.2f %f %f %f %f %f %f %f %f %s %s %s %s\n' % (pdb_id,_chain, qmean_scores.protein_size,fraction_modelled, qmean_scores.energy_Cbeta, qmean_scores.energy_all_atom, qmean_scores.energy_solvation, qmean_scores.energy_torsion,qmean_scores.energy_Cbeta_norm,qmean_scores.energy_all_atom_norm,qmean_scores.energy_solvation_norm,qmean_scores.energy_torsion_norm , qmean_scores.SSE_agree, qmean_scores.ACC_agree, qmean_scores.QMEAN6_score*fraction_modelled,qmean_scores.QMEANnorm6_score*fraction_modelled))

    else:
      out_global.write('%s %s %i %f %f %f %f %f %f %f %f %f %f %s %s %f %f %s %s %s %f %s\n' % (pdb_id,_chain, qmean_scores.protein_size,qmean_scores.fraction_modelled,qmean_scores.fraction_residues_with_sidechains, qmean_scores.energy_Cbeta, qmean_scores.energy_all_atom, qmean_scores.energy_solvation, qmean_scores.energy_torsion,qmean_scores.energy_Cbeta_norm,qmean_scores.energy_all_atom_norm,qmean_scores.energy_solvation_norm,qmean_scores.energy_torsion_norm , qmean_scores.SSE_agree, qmean_scores.ACC_agree, qmean_scores.QMEAN_score*fraction_modelled, qmean_scores.QMEANnorm_score*fraction_modelled,qmean_scores.QMEAN6_score*fraction_modelled,qmean_scores.QMEANnorm6_score*fraction_modelled,qmean_scores.QMEANnorm6_score_trained_on_SwissModel, qmean_scores.QMEANnorm_size_dependent_score,qmean_scores.QMEANnorm6_size_dependent_score))
      if opts.silent_mode == False:
        sys.stdout.write('PDB chain size fraction_modelled fraction_residues_with_sidechains C_beta_pairwise all_atom_pairwise solvation combined_torsion C_beta_pairwise_norm all_atom_pairwise_norm solvation_norm combined_torsion_norm SSE_agree ACC_agree QMEAN4 QMEANnorm QMEAN6_score QMEANnorm6_score QMEANnorm6_score_trained_on_SwissModel QMEANnorm_size_dependent QMEANnorm6_size_dependent\n')
        sys.stdout.write('%s %s %i %f %f %f %f %f %f %f %f %f %f %s %s %f %f %s %s %s %f %s\n' % (pdb_id,_chain, qmean_scores.protein_size,qmean_scores.fraction_modelled,qmean_scores.fraction_residues_with_sidechains, qmean_scores.energy_Cbeta, qmean_scores.energy_all_atom, qmean_scores.energy_solvation, qmean_scores.energy_torsion,qmean_scores.energy_Cbeta_norm,qmean_scores.energy_all_atom_norm,qmean_scores.energy_solvation_norm,qmean_scores.energy_torsion_norm , qmean_scores.SSE_agree, qmean_scores.ACC_agree, qmean_scores.QMEAN_score*fraction_modelled, qmean_scores.QMEANnorm_score*fraction_modelled,qmean_scores.QMEAN6_score*fraction_modelled,qmean_scores.QMEANnorm6_score*fraction_modelled,qmean_scores.QMEANnorm6_score_trained_on_SwissModel, qmean_scores.QMEANnorm_size_dependent_score,qmean_scores.QMEANnorm6_size_dependent_score))

  out_global.close()
  
  #remove temporary pdb (extracted from *.gz)
  if clean_up_flag == 1:
    os.system("rm "+path)
  if opts.silent_mode == False:
    print "QMEAN scores written to:",out_file
  if opts.merge_output == True and opts.analyse_entire_complex == False:
    #calculate Z-score and generate plot:
    Zscore_path=os.path.join(opts.output_dir,pdb_id+"_Zscore.csv")
    plot_path=os.path.join(opts.output_dir,pdb_id+"_plot.png")
    _calculate_Zscore_all(out_file, Zscore_path, plot_path, pdb_id)



def _PdbFileList():
  global qmean_scores

  basename=str(os.getpid())

  structures_with_errors_path=os.path.join(opts.output_dir, "structures_with_errors.txt")
  structures_with_errors_ifile=open(structures_with_errors_path,'a')
  global processed_structures,failed_structures
  processed_structures=[]
  failed_structures=[]

  if opts.identifier !="Txxxx":
    basename =opts.identifier 
  
  if opts.pdb_list != None:  # for list of files from a specific directory (Flo)
    #basename=(os.path.basename(opts.pdb_list)).split('.')[0]
    file_list=open(opts.pdb_list).readlines()
    file_list=[i.split()[0] for i in file_list]

  elif opts.pdb_list == None and opts.pdb_dir != "":  # for all files from a specific directory (e.g. PISCES selections)
    file_list=os.listdir(opts.pdb_dir)

  if opts.calculate_local_scores==True and (opts.psipred==None and opts.sspro==None and opts.predicted_features_directory==None):
    print "Warning: Local score calculation unreliable without PSIPRED and SSpro information"

  out_file=os.path.join(opts.output_dir,"QMEAN_scores_"+basename+".csv")
  out_global=open(out_file, 'w')
  if opts.verbose_output_mode == False:
    out_global.write('PDB chain size fraction_modelled C_beta_pairwise all_atom_pairwise solvation combined_torsion C_beta_pairwise_norm all_atom_pairwise_norm solvation_norm combined_torsion_norm SSE_agree ACC_agree QMEAN_score QMEANnorm_score\n')
  else:
    out_global.write('PDB chain size fraction_modelled fraction_residues_with_sidechains C_beta_pairwise all_atom_pairwise solvation combined_torsion C_beta_pairwise_norm all_atom_pairwise_norm solvation_norm combined_torsion_norm SSE_agree ACC_agree QMEAN4 QMEANnorm QMEAN6_score QMEANnorm6_score QMEANnorm6_score_trained_on_SwissModel QMEANnorm_size_dependent QMEANnorm6_size_dependent\n')
  if opts.silent_mode == False:
    if opts.verbose_output_mode == False:
      sys.stdout.write('PDB chain size fraction_modelled C_beta_pairwise all_atom_pairwise solvation combined_torsion C_beta_pairwise_norm all_atom_pairwise_norm solvation_norm combined_torsion_norm SSE_agree ACC_agree QMEAN_score QMEANnorm_score\n')
    else:
      sys.stdout.write('PDB chain size fraction_modelled fraction_residues_with_sidechains C_beta_pairwise all_atom_pairwise solvation combined_torsion C_beta_pairwise_norm all_atom_pairwise_norm solvation_norm combined_torsion_norm SSE_agree ACC_agree QMEAN4 QMEANnorm QMEAN6_score QMEANnorm6_score QMEANnorm6_score_trained_on_SwissModel QMEANnorm_size_dependent QMEANnorm6_size_dependent\n')

  path=""
  for p in file_list:
    pdb_id=p
    try:
      if len(p)==4 and opts.pdb_list != None:
        PDB_repo=_initialiseRepository(opts.repository)
        ent=PDB_repo.Load(p)  
      elif len(p)==5 and opts.pdb_list != None:
        PDB_repo=_initialiseRepository(opts.repository)
        ent=PDB_repo.Load(p[:4], p[4])
        pdb_id=p[:4].lower()+p[4]
      else:
        path=os.path.join(opts.pdb_dir,p)
        ent=io.LoadPDB(str(path))
        pdb_id=os.path.basename(path)
    except:
      err="Warning: PDB %s not valide" %(p)
      if opts.CASP_output==True:
        os.write(CASP_outstream, p+" X\n")
      failed_structures.append(p)
      continue
      #raise RuntimeError(err)


    if ent.IsValid() == 0 or p == "list" or p == ".svn":
      err="Warning: PDB %s not valide" %(p)
      if opts.CASP_output==True:
        os.write(CASP_outstream, p+" X\n")
      failed_structures.append(p)
      continue
      #raise RuntimeError(err)

    # for set of structures to be assessed either one pair of SSpro/PSIPRED files can be provided or multiple.
    # In the latter case, a directory instead of a file is provided and the path is generated using
    # the directory and the basename (provide to CalculateAgreementTerms() )
    psipred_path=opts.psipred
    sspro_path=opts.sspro
    if opts.psipred!=None and opts.sspro!=None:
    
      is_directory_flag=False
      if opts.pdb_path==None:
        if os.path.isdir(opts.psipred) and os.path.isdir(opts.sspro):
          #print "Directories including predicted sequence featured have been provided."
          # assumption: the basename for the files is provided and storsed in pdb_id
          is_directory_flag=True

      if is_directory_flag==True:
        psipred_path=os.path.join(opts.psipred,pdb_id+".horiz")
        sspro_path=os.path.join(opts.sspro,pdb_id+".sspro_acc")


    if opts.use_occupancy_flag == True:
      ent=ent.Select("occ!=0")

    #add Cbetas if missing:
    qmean.ConstructCBetas(ent, True)

    try:
      if len(ent.chains) > 1 and opts.chain=='' and opts.use_only_first_chain==False and opts.analyse_entire_complex==False:
        if opts.silent_mode == False:
          print "%s chains found for %s..." %(len(ent.chains),pdb_id)
        for c in ent.chains:
          #ent=io.LoadPDB(str(path), c.name)
          SEL='peptide=true and cname='+c.name
          ent_c=ent.Select(SEL)
          #print c.name, ent_c.GetResidueCount()
          if ent_c.GetResidueCount() == 0:
            if opts.verbose_output_mode == True:
              print "Warning: chain %s contains no amino acids (DNA?)" %(c.name)
            continue

          if opts.additional_chain_scores==True:
            qmean_scores=global_scores.GlobalScores((ent_c,ent), psipred_path, sspro_path, stat_pot, opts.QMEAN5, dssp_path=opts.dssp_bin)
          else:
            qmean_scores=global_scores.GlobalScores(ent_c, psipred_path, sspro_path, stat_pot, opts.QMEAN5, dssp_path=opts.dssp_bin)

          #calculate Z-score and generate plot:
          Zscore_path=os.path.join(opts.output_dir,pdb_id+c.name+"_Zscore.csv")
          plot_path=os.path.join(opts.output_dir,pdb_id+c.name+"_plot.png")
          if opts.merge_output == False:
            _calculate_Zscore_single(Zscore_path, plot_path, pdb_id+c.name)

          if opts.calculate_local_scores==True:
            predicted_features_flag=1
            if opts.psipred==None and opts.sspro==None and opts.predicted_features_directory==None:
              predicted_features_flag=0
            try:
              QMEANlocal=local_scores.LocalScores(ent_c, predicted_features_flag, opts.identifier, p, opts.window_size, opts.normalisation_flag, opts.spherical_smoothing_mode, True, opts.QMEAN5)
              if opts.CASP_output==True:
                local_score_string=QMEANlocal.GetCaspVector(qmean_scores.sequence_length)
              #local scores output
              path=os.path.join(opts.output_dir,"energy_profiles/"+pdb_id+c.name+"_local_energy_profile_all_terms.csv")
              if os.path.exists(path):
                os.remove(path)
              local_out_file=os.open(path, os.O_CREAT|os.O_RDWR)
              QMEANlocal.WriteOutput(local_out_file,1)  #0:S-score, 1:RMSD
              if opts.save_Bfactor_colored_PDB==True:
                QMEANlocal.ColorPDB(path_PDB)
              os.close(local_out_file)
              if opts.plot_energy_profile==True:
                _energy_profile(path,os.path.join(opts.output_dir,"energy_profile_plots/"+pdb_id+c.name+"_local_energy_profile_QMEANlocal.png"), pdb_id+c.name)
            except Exception, e:
              raise RuntimeError("Error in local score calculation:", e)

          _chain=c.name
          if _chain=='' or _chain==" ":
              _chain="_"

          fraction_modelled=1.0
          if opts.fraction_modelled_flag == True:
            fraction_modelled=qmean_scores.fraction_modelled
            out_global.write('%s %s %i %.2f %f %f %f %f %f %f %f %f %s %s %s %s\n' % (pdb_id,_chain,qmean_scores.protein_size,fraction_modelled, qmean_scores.energy_Cbeta, qmean_scores.energy_all_atom, qmean_scores.energy_solvation, qmean_scores.energy_torsion,qmean_scores.energy_Cbeta_norm,qmean_scores.energy_all_atom_norm,qmean_scores.energy_solvation_norm,qmean_scores.energy_torsion_norm , qmean_scores.SSE_agree, qmean_scores.ACC_agree, qmean_scores.QMEAN6_score*fraction_modelled,qmean_scores.QMEANnorm6_score*fraction_modelled))
            if opts.silent_mode == False:
              sys.stdout.write('%s %s %i %.2f %f %f %f %f %f %f %f %f %s %s %s %s\n' % (pdb_id,_chain,qmean_scores.protein_size,fraction_modelled, qmean_scores.energy_Cbeta, qmean_scores.energy_all_atom, qmean_scores.energy_solvation, qmean_scores.energy_torsion,qmean_scores.energy_Cbeta_norm,qmean_scores.energy_all_atom_norm,qmean_scores.energy_solvation_norm,qmean_scores.energy_torsion_norm , qmean_scores.SSE_agree, qmean_scores.ACC_agree, qmean_scores.QMEAN6_score*fraction_modelled,qmean_scores.QMEANnorm6_score*fraction_modelled))
          else:
            out_global.write('%s %s %i %f %f %f %f %f %f %f %f %f %f %s %s %f %f %s %s %s %f %s\n' % (pdb_id,_chain, qmean_scores.protein_size,fraction_modelled,qmean_scores.fraction_residues_with_sidechains, qmean_scores.energy_Cbeta, qmean_scores.energy_all_atom, qmean_scores.energy_solvation, qmean_scores.energy_torsion,qmean_scores.energy_Cbeta_norm,qmean_scores.energy_all_atom_norm,qmean_scores.energy_solvation_norm,qmean_scores.energy_torsion_norm , qmean_scores.SSE_agree, qmean_scores.ACC_agree, qmean_scores.QMEAN_score*fraction_modelled, qmean_scores.QMEANnorm_score*fraction_modelled,qmean_scores.QMEAN6_score*fraction_modelled,qmean_scores.QMEANnorm6_score*fraction_modelled,qmean_scores.QMEANnorm6_score_trained_on_SwissModel, qmean_scores.QMEANnorm_size_dependent_score,qmean_scores.QMEANnorm6_size_dependent_score))
            if opts.silent_mode == False:
              sys.stdout.write('%s %s %i %f %f %f %f %f %f %f %f %f %f %s %s %f %f %s %s %s %f %s\n' % (pdb_id,_chain, qmean_scores.protein_size,qmean_scores.fraction_modelled,qmean_scores.fraction_residues_with_sidechains, qmean_scores.energy_Cbeta, qmean_scores.energy_all_atom, qmean_scores.energy_solvation, qmean_scores.energy_torsion,qmean_scores.energy_Cbeta_norm,qmean_scores.energy_all_atom_norm,qmean_scores.energy_solvation_norm,qmean_scores.energy_torsion_norm , qmean_scores.SSE_agree, qmean_scores.ACC_agree, qmean_scores.QMEAN_score, qmean_scores.QMEANnorm_score,qmean_scores.QMEAN6_score,qmean_scores.QMEANnorm6_score,qmean_scores.QMEANnorm6_score_trained_on_SwissModel, qmean_scores.QMEANnorm_size_dependent_score,qmean_scores.QMEANnorm6_size_dependent_score))


      #FOR CASP you are here
      else:
        #add Cbetas if missing:
        qmean.ConstructCBetas(ent, True)
        l=list(ent.chains)
        _chain=" "
        #check if chain has been specified:
        if opts.chain != '':
          _chain=opts.chain
        else: #take first chain
          for i in range(len(l)):
            _chain=((l[i]).name)
            ent_c=ent.Select("peptide=true and cname=\""+_chain+"\"")
            if ent_c.GetResidueCount() != 0:
              break

        ent=ent.Select("peptide=true")
        full_ent=ent
        if ent.GetResidueCount() == 0:
          print "Warning: %s contains no amino acids (DNA?)" %(pdb_id)
          continue
        if opts.analyse_entire_complex==False:
          #ent=io.LoadPDB(str(path), _chain)  
          ent=ent.Select("peptide=true and cname=\""+_chain+"\"")
          qmean_scores=global_scores.GlobalScores(ent, psipred_path, sspro_path, stat_pot, opts.QMEAN5, dssp_path=opts.dssp_bin)

        else:
          if opts.additional_chain_scores==True:
            if (opts.psipred!=None and opts.sspro!=None):
              qmean_scores=global_scores.GlobalScores((ent,full_ent), opts.psipred, opts.sspro, stat_pot, opts.QMEAN5, dssp_path=opts.dssp_bin)
            elif opts.predicted_features_directory!=None:
              qmean_scores=global_scores.GlobalScores((ent,full_ent), opts.predicted_features_directory, opts.predicted_features_directory, stat_pot, opts.QMEAN5, dssp_path=opts.dssp_bin)
            else:
              qmean_scores=global_scores.GlobalScores((ent,full_ent),None,None,stat_pot, dssp_path=opts.dssp_bin)
          else:
            if (opts.psipred!=None and opts.sspro!=None):
              qmean_scores=global_scores.GlobalScores(ent, opts.psipred, opts.sspro, stat_pot, opts.QMEAN5, dssp_path=opts.dssp_bin)
            elif opts.predicted_features_directory!=None:
              qmean_scores=global_scores.GlobalScores(ent, opts.predicted_features_directory, opts.predicted_features_directory, stat_pot, opts.QMEAN5, dssp_path=opts.dssp_bin)
            else:
              qmean_scores=global_scores.GlobalScores(ent,None,None,stat_pot, dssp_path=opts.dssp_bin)
            #raise RuntimeError("For multichain proteins, please specify a directory containing PSIPRED and SSpro prictions in the format MD5.horiz and MD5.sspro_acc")


        #calculate Z-score and generate plot:
        Zscore_path=os.path.join(opts.output_dir,pdb_id+"_Zscore.csv")
        plot_path=os.path.join(opts.output_dir,pdb_id+"_plot.png")
        _calculate_Zscore_single(Zscore_path, plot_path, pdb_id)


        fraction_modelled=1.0
        if opts.fraction_modelled_flag == True:
          fraction_modelled=qmean_scores.fraction_modelled
        if opts.calculate_local_scores!=True:
          if opts.CASP_output==True:
            os.write(CASP_outstream, "%s %.3f\n" %(p, qmean_scores.QMEANnorm6_score*fraction_modelled))
        if opts.calculate_local_scores==True:
            predicted_features_flag=1
            if opts.psipred==None and opts.sspro==None and opts.predicted_features_directory==None:
              print "Warning: Local score calculation unreliable without PSIPRED and SSpro information"
              predicted_features_flag=0
            try:
              QMEANlocal=local_scores.LocalScores(ent, predicted_features_flag, opts.identifier, p, opts.window_size, opts.normalisation_flag, opts.spherical_smoothing_mode, True, opts.QMEAN5)
              if opts.CASP_output==True:
                local_score_string=QMEANlocal.GetCaspVector(qmean_scores.sequence_length)
                fixed_width_string=format_to_fixed_width_string("%s %.3f %s\n" %(p, qmean_scores.QMEANnorm6_score*fraction_modelled, local_score_string))
                os.write(CASP_outstream, fixed_width_string)

              #local scores output
              path=""
              model_name=""
              if  opts.use_only_first_chain==True or opts.analyse_entire_complex==True:
                model_name=pdb_id
              else:
                model_name=pdb_id+"_"+_chain
              #chain_f=open(os.path.join(opts.output_dir,"chain_list.txt"),"w")
              if opts.analyse_entire_complex== True:
                chain_list=[]
                file_list=[]
                for c in ent.chains:
                  model_name=pdb_id+"_"+c.name
                  if c.name==" " or len(ent.chains)==1:
                    model_name=pdb_id
                  chain_list.append(c.name)
                  #chain_f.write(c.name+";")
                  path=os.path.join(opts.output_dir,"energy_profiles/"+model_name+"_local_energy_profile_all_terms.csv")
                  file_list.append(path)

                path_PDB=os.path.join(opts.output_dir,"colored_PDB/"+model_name+"_QMEANlocal_in_Bfactor_column.pdb")
                QMEANlocal.WriteOutputComplexes(chain_list,file_list,1)  #0:S-score, 1:RMSD
                if opts.save_Bfactor_colored_PDB==True:
                  QMEANlocal.ColorPDB(path_PDB)

                for i in range(len(chain_list)):
                  model_name=pdb_id+"_"+chain_list[i]
                  if chain_list[i]==" " or len(len(chain_list))==1:
                    model_name=pdb_id
                  if opts.plot_energy_profile==True:
                    _energy_profile(file_list[i],os.path.join(opts.output_dir,"energy_profile_plots/"+model_name+"_local_energy_profile_QMEANlocal.png"), model_name)

              else:
                path=os.path.join(opts.output_dir,"energy_profiles/"+pdb_id+"_local_energy_profile_all_terms.csv")
                path_PDB=os.path.join(opts.output_dir,"colored_PDB/"+pdb_id+"_QMEANlocal_in_Bfactor_column.pdb")
                if os.path.exists(path):
                  os.remove(path)
                local_out_file=os.open(path, os.O_CREAT|os.O_RDWR)
                QMEANlocal.WriteOutput(local_out_file,1)  #0:S-score, 1:RMSD
                if opts.save_Bfactor_colored_PDB==True:
                  QMEANlocal.ColorPDB(path_PDB)
                os.close(local_out_file)
                if opts.plot_energy_profile==True:
                  _energy_profile(path,os.path.join(opts.output_dir,"energy_profile_plots/"+pdb_id+"_local_energy_profile_QMEANlocal.png"), pdb_id)

            except Exception, e:
              local_score_string=qmean_scores.sequence_length*"X "
              fixed_width_string=format_to_fixed_width_string("%s %.3f %s\n" %(p, qmean_scores.QMEANnorm6_score, local_score_string))
              if opts.CASP_output==True:
                os.write(CASP_outstream, fixed_width_string)
              raise RuntimeError("Error in local score calculation:", e)

        # for complexes -> set chain name to # chains
        if opts.analyse_entire_complex==True:
          _chain=str(ent.GetChainCount())
        if _chain=='':
          _chain="_"

        if opts.verbose_output_mode == False:
          out_global.write('%s %s %i %.2f %f %f %f %f %f %f %f %f %s %s %s %s\n' % (pdb_id,_chain,qmean_scores.protein_size,fraction_modelled, qmean_scores.energy_Cbeta, qmean_scores.energy_all_atom, qmean_scores.energy_solvation, qmean_scores.energy_torsion,qmean_scores.energy_Cbeta_norm,qmean_scores.energy_all_atom_norm,qmean_scores.energy_solvation_norm,qmean_scores.energy_torsion_norm , qmean_scores.SSE_agree, qmean_scores.ACC_agree, qmean_scores.QMEAN6_score*fraction_modelled,qmean_scores.QMEANnorm6_score*fraction_modelled))
          if opts.silent_mode == False:
            sys.stdout.write('%s %s %i %.2f %f %f %f %f %f %f %f %f %s %s %s %s\n' % (pdb_id,_chain,qmean_scores.protein_size,fraction_modelled, qmean_scores.energy_Cbeta, qmean_scores.energy_all_atom, qmean_scores.energy_solvation, qmean_scores.energy_torsion,qmean_scores.energy_Cbeta_norm,qmean_scores.energy_all_atom_norm,qmean_scores.energy_solvation_norm,qmean_scores.energy_torsion_norm , qmean_scores.SSE_agree, qmean_scores.ACC_agree, qmean_scores.QMEAN6_score*fraction_modelled,qmean_scores.QMEANnorm6_score*fraction_modelled))
        else:
          out_global.write('%s %s %i %f %f %f %f %f %f %f %f %f %f %s %s %f %f %s %s %s %f %s\n' % (pdb_id,_chain, qmean_scores.protein_size,fraction_modelled,qmean_scores.fraction_residues_with_sidechains, qmean_scores.energy_Cbeta, qmean_scores.energy_all_atom, qmean_scores.energy_solvation, qmean_scores.energy_torsion,qmean_scores.energy_Cbeta_norm,qmean_scores.energy_all_atom_norm,qmean_scores.energy_solvation_norm,qmean_scores.energy_torsion_norm , qmean_scores.SSE_agree, qmean_scores.ACC_agree, qmean_scores.QMEAN_score*fraction_modelled, qmean_scores.QMEANnorm_score*fraction_modelled,qmean_scores.QMEAN6_score*fraction_modelled,qmean_scores.QMEANnorm6_score*fraction_modelled,qmean_scores.QMEANnorm6_score_trained_on_SwissModel, qmean_scores.QMEANnorm_size_dependent_score,qmean_scores.QMEANnorm6_size_dependent_score))
          if opts.silent_mode == False:
            sys.stdout.write('%s %s %i %f %f %f %f %f %f %f %f %f %f %s %s %f %f %s %s %s %f %s\n' % (pdb_id,_chain, qmean_scores.protein_size,qmean_scores.fraction_modelled,qmean_scores.fraction_residues_with_sidechains, qmean_scores.energy_Cbeta, qmean_scores.energy_all_atom, qmean_scores.energy_solvation, qmean_scores.energy_torsion,qmean_scores.energy_Cbeta_norm,qmean_scores.energy_all_atom_norm,qmean_scores.energy_solvation_norm,qmean_scores.energy_torsion_norm , qmean_scores.SSE_agree, qmean_scores.ACC_agree, qmean_scores.QMEAN_score, qmean_scores.QMEANnorm_score,qmean_scores.QMEAN6_score,qmean_scores.QMEANnorm6_score,qmean_scores.QMEANnorm6_score_trained_on_SwissModel, qmean_scores.QMEANnorm_size_dependent_score,qmean_scores.QMEANnorm6_size_dependent_score))
  
    except Exception, e:
      failed_structures.append(p)
      structures_with_errors_ifile.write('%s\t%s\n' %(p,e))
      if opts.silent_mode == False:
        print e

  out_global.close()
  if opts.silent_mode == False:
    print "QMEAN scores written to:",out_file
  if opts.merge_output == True:
    #calculate Z-score and generate plot:
    Zscore_path=os.path.join(opts.output_dir,basename+"_Zscore.csv")
    plot_path=os.path.join(opts.output_dir,basename+"_plot.png")
    _calculate_Zscore_all(out_file, Zscore_path, plot_path, basename)



def _calculate_Zscore_single(Zscore_path, plot_path, pdb_id):
  if opts.no_zscores==True:
    return
  R_abs_path=settings.Locate('R', env_name='R_EXECUTABLE')
  if (opts.psipred!=None and opts.sspro!=None) or opts.predicted_features_directory!=None:
    os.system(R_abs_path+" < "+conf.QMEAN_ZSCORE_DATA_DIR+"/calculate_Zscore_from_QMEANnorm_scores.R %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s --no-save  >/dev/null" %(qmean_scores.QMEANnorm6_score, qmean_scores.energy_Cbeta_norm,qmean_scores.energy_all_atom_norm,qmean_scores.energy_solvation_norm,qmean_scores.energy_torsion_norm, qmean_scores.SSE_agree, qmean_scores.ACC_agree, qmean_scores.protein_size, Zscore_path, plot_path, pdb_id, opts.Membrane_potentials, opts.save_thumbnail, opts.silent_mode, opts.no_plots, opts.label_structures, opts.analyse_entire_complex, conf.QMEAN_ZSCORE_DATA_DIR, opts.quality_R_plots, opts.slider_flag))
  else:
    #print conf.QMEAN_ZSCORE_DATA_DIR+"/calculate_Zscore_from_QMEANnorm_scores.R"
    os.system(R_abs_path+" < "+conf.QMEAN_ZSCORE_DATA_DIR+"/calculate_Zscore_from_QMEANnorm_scores.R %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s --no-save >/dev/null" %(qmean_scores.QMEANnorm_score, qmean_scores.energy_Cbeta_norm,qmean_scores.energy_all_atom_norm,qmean_scores.energy_solvation_norm,qmean_scores.energy_torsion_norm, qmean_scores.SSE_agree, qmean_scores.ACC_agree, qmean_scores.protein_size, Zscore_path, plot_path, pdb_id, opts.Membrane_potentials, opts.save_thumbnail, opts.silent_mode, opts.no_plots, opts.label_structures, opts.analyse_entire_complex, conf.QMEAN_ZSCORE_DATA_DIR, opts.quality_R_plots, opts.slider_flag))


def _calculate_Zscore_all(QMEAN_scores_file, Zscore_path, plot_path, pdb_id):
  if opts.no_zscores==True:
    return
  R_abs_path=settings.Locate('R', env_name='R_EXECUTABLE')
  QMEAN_mode="QMEANnorm"
  if (opts.psipred!=None and opts.sspro!=None) or opts.predicted_features_directory!=None:
    QMEAN_mode="QMEANnorm6"
  os.system(R_abs_path+" < "+conf.QMEAN_ZSCORE_DATA_DIR+"/calculate_Zscore_from_QMEANnorm_scores_table.R %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s --no-save >/dev/null" %(QMEAN_scores_file, QMEAN_mode, Zscore_path, plot_path, pdb_id, opts.Membrane_potentials, opts.save_thumbnail, opts.silent_mode, opts.no_plots, opts.assign_experimental_method, opts.label_structures, opts.analyse_entire_complex, conf.QMEAN_ZSCORE_DATA_DIR, opts.quality_R_plots, opts.slider_flag))


def _energy_profile(results_table, out_filename, model_name):
  R_abs_path=settings.Locate('R', env_name='R_EXECUTABLE')
  os.system(R_abs_path+" < "+conf.QMEAN_ZSCORE_DATA_DIR+"/draw_energy_profile.R %s %s %s %s %s %s --no-save >/dev/null" %(results_table, out_filename, model_name, opts.quality_R_plots,opts.local_distance_measure, opts.save_thumbnail))  #0:S-score, 1:RMSD




def _CASP_header():
  qmode=1
  if opts.calculate_local_scores==True:
    qmode=2
  header=""
  if opts.method=="QMEANclust":
    header="""PFRMAT QA
TARGET %s
AUTHOR %s
METHOD QMEANclust combines structural desity information provided by
METHOD the ensemble of models with the QMEAN scoring function.
METHOD Clustering is performed based on the models pairwise GDT_TS score.
METHOD Only a fraction of the models with the best QMEAN score is
METHOD used in order to derive the structural desity information.
METHOD Reference:
METHOD Benkert, P., Schwede, T. and Tosatto, S.C.E. (2009). 
METHOD \"QMEANclust: Estimation of protein model quality by "
METHOD combining a composite scoring function with structural 
METHOD density information."\"; BMC Struct Biol. 2009 May 20;9:35.
MODEL 1
QMODE %s
""" %(opts.identifier,opts.method,qmode)

  elif opts.method=="QMEANdist":
    header="""PFRMAT QA
TARGET %s
AUTHOR %s
METHOD QMEAN combined with distance constrains from templates
METHOD Reference:
METHOD Benkert, P., Tosatto, S.C.E. and Schomburg, D. (2008). 
METHOD \"QMEAN: A comprehensive scoring function for model 
METHOD quality assessment.\"; Proteins, 71(1):261-277.
MODEL 1
QMODE %s
""" %(opts.identifier,opts.method,qmode)

  else:
    header="""PFRMAT QA
TARGET %s
AUTHOR %s
METHOD QMEAN is a composite scoring function consisting of six
METHOD structural descriptors: The local geometry is analysed by a
METHOD torsion angle potential over three consecutive amino acids.
METHOD A distance-dependent pairwise Cbeta potential as well as
METHOD an all-atom potential with 167 atom types are used to
METHOD assess long-range interactions. A solvation potential describes
METHOD the burial status of the residues. Two simple terms describing
METHOD the agreement of predicted and calculated secondary structure
METHOD and solvent accessibility, respectively, are also included.
METHOD Reference:
METHOD Benkert, P., Tosatto, S.C.E. and Schomburg, D. (2008). 
METHOD \"QMEAN: A comprehensive scoring function for model 
METHOD quality assessment.\"; Proteins, 71(1):261-277.
MODEL 1
QMODE %s
""" %(opts.identifier,opts.method,qmode)

  return header



def format_to_fixed_width_string(s,width=80):
  new_s=""
  pos=width
  while(len(s) >=width):
    while(s[pos]!=" " and pos!=len(s)-1):
      pos+=1
    new_s=new_s+s[:pos]+"\n"
    s=s[pos:]
    pos=width
  new_s+=s
  return new_s



#MAIN THIS IS (Yoda):
#try:
opts=_ParseOptions()

if opts.Membrane_potentials==True:
  stat_pot=statistical_potentials.StatisticalPotentials(conf.MEMBRANE_PROTEIN_POTENTIALS_DIR)
  if opts.silent_mode == False:
    print "Membrane-derived statistical potentials used"
else:
  stat_pot=statistical_potentials.StatisticalPotentials(conf.POTENTIALS_DIR)
  #stat_pot=statistical_potentials.StatisticalPotentials("/import/bc2/home/schwede/benkert/qmean/data/qmean/potentials/repaired_Cbeta")

#CASP output:
if opts.CASP_output==True:
  global CASP_outstream
  out_filename=os.path.join(opts.output_dir,"QMEAN_CASP_%s.txt" %opts.identifier)
  
  if os.path.exists(out_filename):
    os.remove(out_filename)
  CASP_outstream=os.open(out_filename, os.O_CREAT|os.O_WRONLY)
  os.write(CASP_outstream, _CASP_header())

if opts.calculate_local_scores==True:
  os.system("mkdir -p "+os.path.join(opts.output_dir,"energy_profiles"))
if opts.plot_energy_profile==True:
  os.system("mkdir -p "+os.path.join(opts.output_dir,"energy_profile_plots"))
if opts.save_Bfactor_colored_PDB==True:
  os.system("mkdir -p "+os.path.join(opts.output_dir,"colored_PDB"))

# depending on option either call function for handling single PDB or multiple:
if opts.pdb_path == None and (opts.pdb_list != None or opts.pdb_dir != ""): # for list of PDB identifiers (e.g. PISCES selections)
  _PdbFileList()

elif opts.pdb_path != None:
  if len(opts.pdb_path) == 4:
    pdb_id=(opts.pdb_path).lower()
    chain=opts.chain
  elif len(opts.pdb_path) == 5:
    pdb_id=(opts.pdb_path)[:-1].lower()
    chain=(opts.pdb_path)[-1].upper()
  else:
    pdb_id=opts.pdb_path
    chain=opts.chain
  _SinglePdb(pdb_id, chain)

else:
  print "Either specify a single structure or multiple structures (list or directory)"
  
if opts.CASP_output==True:
  os.write(CASP_outstream, "END")
  os.close(CASP_outstream)

#workaround for adding global per-chain scores:
if opts.additional_chain_scores==True and opts.analyse_entire_complex == True:
  opts.analyse_entire_complex = False
  opts.calculate_local_scores=False
  opts.merge_output = True
  opts.output_dir+="/chains"
  os.system("mkdir -p "+opts.output_dir)
  _SinglePdb(pdb_id, chain)

#except Exception, e:
#  raise e  #print e
#










