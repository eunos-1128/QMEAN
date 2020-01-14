import os
import random
from ost.table import Table
from qmean import LocalScorer

#we need to have some training data
training_tab = Table()

#The table must have a field called "target" => that's the target value
#we're training for.
training_tab.AddCol("target",'f')

#The table also must have a field called AA to train the AA specific linear
#combinations
training_tab.AddCol('AA','s')

#Let's add some amazing features, that should predict the target
features = ["feature_1", "feature_2", "feature_3"]
training_tab.AddCol(features[0],'f')
training_tab.AddCol(features[1],'f')
training_tab.AddCol(features[2],'f')

one_letter_codes = ['A', 'R', 'N', 'D', 'Q', 'E', 'K', 'S', 'C', 'M', 
                    'W', 'Y', 'T', 'V', 'I', 'L', 'G', 'P', 'H', 'F']

#let's add some random training data for every amino acid

for olc in one_letter_codes:
  
  for i in range(100):
    
    row = dict()
    row["AA"] = olc
    row["target"] = random.random()
    row[features[0]] = random.random()
    row[features[1]] = random.random()
    row[features[2]] = random.random()
    training_tab.AddRow(row)

#Let's setup a LocalScorer
local_scorer = LocalScorer()

#and train it with the random data
local_scorer.TrainNewType(features,training_tab,"random")

data_one = dict()
data_two = dict()

#let's get some local score given the data dicts
data_one[features[0]] = 1.0
data_one[features[1]] = 1.0
data_one[features[2]] = 1.0

data_two[features[0]] = 1.0
data_two[features[1]] = 1.0
data_two[features[2]] = float("NaN")

print "Local score for alanine when all features are defined:"
print local_scorer.GetLocalScore("random",'A',data_one)

print "Local score for alanine when one of the features is NaN:"
print local_scorer.GetLocalScore("random",'A',data_two)

#save down the scorer...
local_scorer.Save("local_scorer.dat")
