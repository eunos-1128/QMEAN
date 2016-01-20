#requires mlabwrap and matlab! download mlabwrap from sourceforge and
#try to run "python setup.py install" If you get errors like strange 
#int conversion stuff, modify setup.py according to http://obasic.net/how-to-install-mlabwrap-on-windows
#If the matlab engine can't be started, make sure that csh is installed!
#
#In case of the local scorer, the ost table MUST contain the column AA (AminoAcid one letter code)



from ost.table import *
import numpy as np
import sys
import traceback

class Scorer():

  def __init__(self):

    self.weights=None
    self.trained_features=None

  def TrainWeights(self, table, feature_names, method='leastsqr'):

    data=list()
    target=list()
    feature_indices=list()

    try:
      target_idx=table.GetColIndex('target')
    except:
      print "training table requires a column named target to train against!"
      sys.exit(0)

    for f in feature_names:
      try:
        feature_indices.append(table.GetColIndex(f))
      except:
        print "the feature\"",f,"\" could not be found as a column in the provided table"
        sys.exit(0)

    for r in table.rows:
      temp=list()
      valid_features=True
      valid_target=True
      for idx in feature_indices:
        value=r[idx]
        if value!=value or value==None:
          valid_features=False
          break
        temp.append(value)

      target_value=r[target_idx]

      if target_value!=target_value or target_value==None:
        valid_target=False

      if valid_features and valid_target:  
        data.append(temp)
        target.append(r[target_idx])
    

    numpy_data=np.zeros((len(data),len(feature_names)))
    numpy_target_data=np.zeros((len(data),1))

    for i in range(len(data)):
      numpy_target_data[i]=target[i]
      for j in range(len(feature_names)):
        numpy_data[i][j]=data[i][j]

    if len(data)<100:
      print 'not enough datapoints! set all weights to 0!!!'
      self.weights=[0]*(len(feature_names)+1)
      return

    if method=='robust':
      try:
       from mlabwrap import mlab
      except ImportError:
        print "If method of feature combination is set to robust, mlabwrap and matlab are necessary! Could not import mlabwrap!"
        raise

      self.weights=mlab.robustfit(numpy_data, numpy_target_data)
      print "trained weights: ",self.weights
      self.trained_features=feature_names

    elif method=='leastsqr':
      try:
        from scipy.linalg import lstsq
      except ImportError:
        print "If method of feature combination is set to leastsqr, scipy.linalg.lstsq is necessary! Could not import it!"

      temp=np.ones((numpy_data.shape[0],numpy_data.shape[1]+1))
      temp[:,1:]=numpy_data
      numpy_data=temp
      self.weights=lstsq(numpy_data, numpy_target_data)[0]
      print "trained weights: ",self.weights
      self.trained_features=feature_names

    else:
      raise ValueError('Only support leastsqr and robust for feature combinations, not '+method)


  def GetScore(self, data):
    
    feature_values=np.ones(len(self.trained_features)+1)
    set_values=0
    
    for i, f in enumerate(self.trained_features):
      try:
        feature_values[i+1]=data[f]
        set_values+=1
      except:
        print f+' is not trained! Return NaN for local score!'
        return float('NaN')

    if set_values!=len(self.trained_features):
      print 'There is one or more features missing! Return NaN for local score!'
      return float('NaN')

    return float(np.dot(feature_values,self.weights))

  def GetWeights(self):
    return self.weights

  def Save(self, stream_or_filename):
    if not hasattr(stream_or_filename, 'write'):
      stream=open(stream_or_filename, 'wb')
    else:
      stream=stream_or_filename
    cPickle.dump(self, stream, cPickle.HIGHEST_PROTOCOL)

  @staticmethod
  def Load(stream_or_filename):
    if not hasattr(stream_or_filename, 'read'):
      stream=open(stream_or_filename, 'rb')
    else:
      stream=stream_or_filename
    return cPickle.load(stream)


class LocalScorer():

  def __init__(self):
    self.scorer=dict()
    self.feature_combinations=dict()
    self.trained_types=list()
    self.AANames=['A', 'R', 'N', 'D', 'Q', 'E', 'K', 'S', 'C', 'M', 'W', 'Y', 'T', 'V', 'I', 'L', 'G', 'P', 'H', 'F']
    self.last_used_features=list()

  def TrainNewType(self, features, tab, residue_type, method='robust'):

    import itertools

    #do some checks
    try:
      tab.GetColIndex('target')
    except:
      raise ValueError('training table requires a column named target to train against!')
    for f in features:
      try:
        tab.GetColIndex(f)
      except:
        raise ValueError('cannot train with feature '+f+', since it is not present in the training table!')
 
    self.feature_combinations[residue_type]=list()

    #build all possible permutations of the features
    for i in range(1,len(features)+1):
      combinations=itertools.combinations(features,i)
      for c in combinations:
        self.feature_combinations[residue_type].append(c)

    #select the tabs for single amino acids
    subtabs=dict()
    for AA in self.AANames:
      subtabs[AA]=tab.Select('AA='+AA)

    self.scorer[residue_type]=dict()

    for AA in self.AANames:
      self.scorer[residue_type][AA]=list()
      for c in self.feature_combinations[residue_type]:
        print "train tab with features ", c," for residue_type ",residue_type," and amino acid ",AA
        self.scorer[residue_type][AA].append(Scorer())
        self.scorer[residue_type][AA][-1].TrainWeights(subtabs[AA], list(c), method)

    self.trained_types.append(residue_type)


  def GetLocalScore(self, residue_type, AA, data):

    if residue_type not in self.trained_types:
      print "provided residue type is not trained! returned NaN for local score!"
      return float('NaN')
    if AA not in self.AANames:
      print "provided amino acid name is not trained! returned NaN for local score!"
      return float('NaN')

    #remove NaN and None values
    processed_data=dict()
    for k,v in data.iteritems():
      if v!=v or v==None:
        continue
      processed_data[k]=v

    #find the right feature combination to get corresponding scorer
    for i, c in enumerate(self.feature_combinations[residue_type]):
      if len(c)!=len(processed_data):
        continue

      match=True

      for k,v in processed_data.iteritems():
        if k not in c:
          match=False

      if match==False:
        continue

      self.last_used_features=list(c)
      return self.scorer[residue_type][AA][i].GetScore(processed_data) 
    print "could not find appropriate scorer for provided data. return NaN for local score"
    return float('NaN')

  def GetWeights(self, residue_type, AA, features):

    if residue_type not in self.trained_types:
      raise ValueError("provided residue type is not trained")
    if AA not in self.AANames:
      raise ValueError("provided amino acid name is not trained")

    for i, c in enumerate(self.feature_combinations[residue_type]):
      if len(c)!=len(features):
        continue

      match=True

      for k in features:
        if k not in c:
          match=False

      if match==False:
        continue

      return self.scorer[residue_type][AA][i].GetWeights() 

  def GetLastUsedFeatures(self):
    return


  def Save(self, stream_or_filename):
    if not hasattr(stream_or_filename, 'write'):
      stream=open(stream_or_filename, 'wb')
    else:
      stream=stream_or_filename
    cPickle.dump(self, stream, cPickle.HIGHEST_PROTOCOL)


  @staticmethod
  def Load(stream_or_filename):
    if not hasattr(stream_or_filename, 'read'):
      stream=open(stream_or_filename, 'rb')
    else:
      stream=stream_or_filename
    return cPickle.load(stream)



class GlobalScorer():

  def __init__(self):
    self.scorer=dict()
    self.feature_combinations=dict()
    self.trained_types=list()

  def TrainNewType(self, features, tab, structure_type, method='robust'):

    import itertools

    #do some checks
    try:
      tab.GetColIndex('target')
    except:
      raise ValueError('training table requires a column named target to train against!')
    for f in features:
      try:
        tab.GetColIndex(f)
      except:
        raise ValueError('cannot train with feature '+f+', since it is not present in the training table!')
 
    self.feature_combinations[structure_type]=list()

    #build all possible permutations of the features
    for i in range(1,len(features)+1):
      combinations=itertools.combinations(features,i)
      for c in combinations:
        self.feature_combinations[structure_type].append(c)

    self.scorer[structure_type]=list()

    for c in self.feature_combinations[structure_type]:
      print "train tab with features ", c," for structure_type ",structure_type
      self.scorer[structure_type].append(Scorer())
      self.scorer[structure_type][-1].TrainWeights(tab, list(c), method)

    self.trained_types.append(structure_type)


  def GetGlobalScore(self, structure_type, data):

    if structure_type not in self.trained_types:
      print "provided structure type is not trained! return NaN for global score!"
      return float('NaN')

    #remove NaN and None values
    processed_data=dict()
    for k,v in data.iteritems():
      if v!=v or v==None:
        continue
      processed_data[k]=v

    #find the right feature combination to get corresponding scorer
    for i, c in enumerate(self.feature_combinations[structure_type]):
      if len(c)!=len(processed_data):
        continue

      match=True

      for k,v in processed_data.iteritems():
        if k not in c:
          match=False

      if match==False:
        continue

      return self.scorer[structure_type][i].GetScore(processed_data) 
    print "could not find appropriate scorer for provided data. return NaN for global score"
    return float('NaN')

  def GetWeights(self, structure_type, features):

    if structure_type not in self.trained_types:
      raise ValueError("provided residue type is not trained")

    for i, c in enumerate(self.feature_combinations[residue_type]):
      if len(c)!=len(features):
        continue

      match=True

      for k in features:
        if k not in c:
          match=False

      if match==False:
        continue

      return self.scorer[structure_type][i].GetWeights() 


  def Save(self, stream_or_filename):
    if not hasattr(stream_or_filename, 'write'):
      stream=open(stream_or_filename, 'wb')
    else:
      stream=stream_or_filename
    cPickle.dump(self, stream, cPickle.HIGHEST_PROTOCOL)


  @staticmethod
  def Load(stream_or_filename):
    if not hasattr(stream_or_filename, 'read'):
      stream=open(stream_or_filename, 'rb')
    else:
      stream=stream_or_filename
    return cPickle.load(stream)

