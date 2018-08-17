import os.path
import ost.settings as bs
import ost
import sys

# solutions with cmake configure_file won't work, as the install directory
# can be adapted when calling make install. This soultion assumes conf.py
# to be in stage/lib64/python2.7/site-packages/qmean
QMEAN_DATA_DIR = os.path.abspath(os.path.join(os.path.dirname( __file__ ), 
                                              '..', '..', '..', '..', 'share'))
if not os.path.exists(QMEAN_DATA_DIR):
  print 'internal configuration error:'
  print '  The QMean data directory "%s" does not exist.' % QMEAN_DATA_DIR
  print '   Make sure QMEAN_DATA_DIR is set correctly'
  sys.exit(-1)
  
POTENTIAL_DIR=os.path.join(QMEAN_DATA_DIR, 'potentials')
REFERENCE_DIR=os.path.join(QMEAN_DATA_DIR, 'reference_values')
SCORER_DIR=os.path.join(QMEAN_DATA_DIR,'scorer')

class SwissmodelSettings:
  def __init__(self):
    self.local_potentials = os.path.join(POTENTIAL_DIR, 'soluble_local_potentials.dat')
    self.global_potentials = os.path.join(POTENTIAL_DIR, 'soluble_global_potentials.dat')
    self.local_scorer = os.path.join(SCORER_DIR, 'local_nn_scorer')
    self.global_scorer = os.path.join(SCORER_DIR, 'promod_global_scorer.dat')
    self.reference_tab = os.path.join(REFERENCE_DIR, 'reference_tab_promod_scorer.txt') 

class MembraneSettings:
  def __init__(self):
    self.local_potentials_soluble = os.path.join(POTENTIAL_DIR, 'soluble_local_potentials.dat')
    self.global_potentials_soluble = os.path.join(POTENTIAL_DIR, 'soluble_global_potentials.dat')
    self.local_potentials_membrane = os.path.join(POTENTIAL_DIR, 'membrane_local_potentials.dat')
    self.global_potentials_membrane = os.path.join(POTENTIAL_DIR, 'membrane_global_potentials.dat')
    self.local_scorer_soluble = os.path.join(SCORER_DIR, 'promod_local_scorer.dat')
    self.local_scorer_membrane = os.path.join(SCORER_DIR, 'membrane_local_scorer.dat')
    self.global_scorer = os.path.join(SCORER_DIR, 'promod_global_scorer.dat')
    self.reference_tab = os.path.join(REFERENCE_DIR, 'reference_tab_promod_scorer.txt')
 
