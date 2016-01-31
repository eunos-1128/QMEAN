#Example code to demonstrate the usage of the container classes

#load required modules
from qmean import *

#define some training targets used later on
training_targets = ["7ODC","2QUD","3T47","4C1A","3ZSJ","4I4T","3H6R","3BSO",
                    "3M0Z","2GPI","3MHP","3V19","2IXS","4G6T","3ZG9","3WSG",
                    "1JH6","2XZE","2OLR","1TIF","3C4S","2CVI","4EET","1LM5",
                    "1MJ5","2FJR","2D3G","3ZNV","3WA2","3WU2","4M7R","2PR7",
                    "3FTJ","1KD8","1HBN","4TRK","2GB4","3HNY","4R7Q","1EAQ",
                    "2O1Q","4DX5","1XAK","5CSM","2XWP","2UWA","3SQZ","4MT8",
                    "3HTU"]


#We define two different statistics in a container and fill them with 
#structural information
s_container = StatisticContainer()
s_container["int"] = InteractionStatistic(0.0,10.0,10,4)
s_container["pack"] = PackingStatistic(5.0,32)

for t in training_targets:
    prot = io.LoadPDB(t,remote=True).Select("peptide=true")
    s_container["int"].Extract(prot,prot)
    s_container["pack"].Extract(prot,prot)

#instead of creating single potentials, we also use the power
#of the potential container
p_container = PotentialContainer()
p_container["int"] = InteractionPotential.Create(s_container["int"])
p_container["pack"] = PackingPotential.Create(s_container["pack"])

#let's get a crambin and print some global energies
crambin = io.LoadPDB("1crn",remote=True).Select("peptide=true")

print "interaction e: ", p_container["int"].GetTotalEnergy(crambin,crambin)
print "packing e: ", p_container["pack"].GetTotalEnergy(crambin,crambin)


#The awesome thing is, that we can dump the full container at once
p_container.Save("potential_container.dat")
#and load it again
new_p_container = PotentialContainer.Load("potential_container.dat")

print "Just loaded the same container again, let's see whether the "
print "energies are consistent"

print "interaction e: ", new_p_container["int"].GetTotalEnergy(crambin,crambin)
print "packing e: ", new_p_container["pack"].GetTotalEnergy(crambin,crambin)
