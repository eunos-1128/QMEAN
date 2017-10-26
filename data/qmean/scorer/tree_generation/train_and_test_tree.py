import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn import metrics
from scipy.stats import pearsonr
import pickle

# Loads the csv files created by tree_data_extraction.py and trains 
# a random forest. It directly gets tested on independent test data 
# using some basic metrics.
# Pearson correlation of prediction vs target and ROC AUC.
# Depending on your version of sklearn, the results might slightly differ
# from the example output below (but should be in the same range!).  


training_file_path = "training_data.csv"
test_file_path = "test_data.csv"
random_forest_path = "disco_tree.pkl"

training_data = open(training_file_path,'r').readlines()

features = [item.strip() for item in training_data[0].strip().split(',')]
required_features = ["qmean4", "avg_max_seqsim", "avg_num_cluster",
                     "avg_variance", "counts", "avg_max_seqid", "local_qmean", 
                     "disco"]

feature_indices = list()
for f in required_features:
  # This will raise an error if f is not there
  feature_indices.append(features.index(f))
target_idx = features.index("lddt")
local_qmean_idx = features.index("local_qmean")
local_disco_idx = features.index("disco")


print "START TO FILL DATA"

# fill it into simple lists of data and kick out any values containing NaN
feature_data = list()
target_data = list()

total_num_data_points = 0
total_num_valid_data_points = 0

for i,line in enumerate(training_data[1:]):

  split_data = line.strip().split(',')

  all_good = True

  temp_data = list()
  for j, f_idx in enumerate(feature_indices):
    val = float(split_data[f_idx])
    if val != val:
      all_good = False
      break
    temp_data.append(val)

  target_value = float(split_data[target_idx])
  if target_value != target_value:
    all_good = False

  if not all_good:
    print "FOUND NaN VALUE IN DATA WITH IDX ", i, "SKIP..."

  else:
    feature_data.append(temp_data)
    target_data.append(target_value)  
    total_num_valid_data_points += 1
  total_num_data_points += 1

print "TOTAL NUM DATA POINTS: ", total_num_data_points
print "TOTAL NUM VALID DATA POINTS: ", total_num_valid_data_points

# dump it into numpy array
feature_array = np.zeros((total_num_valid_data_points, len(required_features))) 
target_array = np.zeros(total_num_valid_data_points)

for i in range(total_num_valid_data_points):
  for j in range(len(required_features)):
    feature_array[(i,j)] = feature_data[i][j]
  target_array[i] = target_data[i]

clf = RandomForestRegressor(n_estimators=10, min_samples_split=100, min_samples_leaf=100,
                            max_features='log2')

print "START TO FIT"
clf.fit(feature_array, target_array)

print "DUMP TO DISK"
pickle.dump(clf,open(random_forest_path,"wb"))

print "START TESTING"

test_data = open(test_file_path,'r').readlines()

qmean_target_values = list()
qmean_values = list()

disco_target_values = list()
disco_values = list()

qmeandisco_target_values = list()
qmeandisco_values = list()

for i,line in enumerate(test_data[1:]):

  split_data = line.strip().split(',')

  all_there = True

  tree_input_data = np.zeros((1, len(feature_indices)))

  for j, f_idx in enumerate(feature_indices):
    val = float(split_data[f_idx])
    if val != val:
      all_there = False
      break
    tree_input_data[(0, j)] = val


  local_qmeandisco_score = float("NaN")
  local_qmean_score = float("NaN")
  local_disco_score = float("NaN")

  local_qmean_score = float(split_data[local_qmean_idx])
  local_disco_score = float(split_data[local_disco_idx])
  local_lddt_score = float(split_data[target_idx])
  
  if all_there:
    # Its all there! lets use the awesome random forest for qmeandisco!
    local_qmeandisco_score = clf.predict(tree_input_data)[0]
  else:
    # QMEANDisCo uses the QMEAN or DisCo value as fallback if any of the 
    # required features is NaN
    if local_qmean_score == local_qmean_score:
      local_qmeandisco_score = local_qmean_score
    elif local_disco_score == local_disco_score:
      local_qmeandisco_score = local_disco_score

  if local_lddt_score == local_lddt_score and local_qmean_score == local_qmean_score:
    qmean_target_values.append(local_lddt_score)
    qmean_values.append(local_qmean_score)
    
  if local_lddt_score == local_lddt_score and local_disco_score == local_disco_score:
    disco_target_values.append(local_lddt_score)
    disco_values.append(local_disco_score)

  if local_lddt_score == local_lddt_score and local_qmeandisco_score == local_qmeandisco_score:
    qmeandisco_target_values.append(local_lddt_score)
    qmeandisco_values.append(local_qmeandisco_score)


def ROCAUC(lddt, score):
  classes = list()
  for item in lddt:
    if item < 0.6:
      classes.append(0)
    else:
      classes.append(1)
  fpr, tpr, thresholds = metrics.roc_curve(classes, score)
  return metrics.auc(fpr, tpr)

print "RANDOM FOREST TEST"
print
print "NUM DATA POINTS:"
print "QMEAN:", len(qmean_values)
print "DISCO:", len(disco_values)
print "QMEANDISCO", len(qmeandisco_values)
print
print "PEARSON CORRELATION"
print "QMEAN:", pearsonr(qmean_values, qmean_target_values)[0]
print "DISCO:", pearsonr(disco_values, disco_target_values)[0]
print "QMEANDISCO:", pearsonr(qmeandisco_values, qmeandisco_target_values)[0]
print
print "ROC AUC:"
print "QMEAN:", ROCAUC(qmean_target_values, qmean_values)
print "DISCO:", ROCAUC(disco_target_values, disco_values)
print "QMEANDISCO:", ROCAUC(qmeandisco_target_values, qmeandisco_values)


# THIS IS THE EXPECTED OUTPUT OF THE ABOVE
# ALL GOOD IF ITS SIMILAR
#
#NUM DATA POINTS:
#QMEAN: 503735
#DISCO: 502515
#QMEANDISCO 503735
#
#PEARSON CORRELATION
#QMEAN: 0.732438837207
#DISCO: 0.824140240048
#QMEANDISCO: 0.865246334716
#
#ROC AUC:
#QMEAN: 0.878284468368
#DISCO: 0.915929981254
#QMEANDISCO: 0.93541767048
