#Example code to demonstrate the usage of setting up an interaction potential
#and get some basic info out of it. You need matplotlib to run this script!

#load required modules
from qmean import *
import matplotlib.pyplot as plt


#define some training targets used later on
training_targets = ["7ODC","2QUD","3T47","4C1A","3ZSJ","4I4T","3H6R","3BSO",
                    "3M0Z","2GPI","3MHP","3V19","2IXS","4G6T","3ZG9","3WSG",
                    "1JH6","2XZE","2OLR","1TIF","3C4S","2CVI","4EET","1LM5",
                    "1MJ5","2FJR","2D3G","3ZNV","3WA2","3WU2","4M7R","2PR7",
                    "3FTJ","1KD8","1HBN","4TRK","2GB4","3HNY","4R7Q","1EAQ",
                    "2O1Q","4DX5","1XAK","5CSM","2XWP","2UWA","3SQZ","4MT8",
                    "3HTU"]


#Create an interaction statistic, considering pairwise
#distances from 0.0 to 10 A with a binsize of 1 A.
#Pairwise interactions from residues separated by less
#than 4 residues are neglected.
stat = InteractionStatistic(0.0,10.0,10,4)


#Let's fill the stat with some data!
for t in training_targets:
    prot = io.LoadPDB(t,remote=True).Select("peptide=true")
    stat.Extract(prot,prot)


#Create a potential object from the extracted statistic
pot = InteractionPotential.Create(stat)

#Save the statistic and potential objects
stat.Save("interaction_statistic.dat")
pot.Save("interaction_pot.dat")


#We're not done yet! Let's play around a bit

#Let's investigate the interaction between a typical hyrogen bond
#donor and acceptor...

pseudo_energies = list()
for i in range(10):
    distance = i
    e = pot.GetEnergy(ChemType.O_SER_G, ChemType.N_HIS_E2, distance)
    pseudo_energies.append(e)

#let's produce an awesome plot and save it down
x = list()
y = list()

for i in range(10):
    x.append(i)
    x.append(i+1)
    y.append(pseudo_energies[i])
    y.append(pseudo_energies[i])

plt.plot(x,y,label="SerO - HisNE2",color='#004586')
plt.plot([0,10],[0,0],'k')
plt.xlabel("distance")
plt.ylabel("pseudo energy")
plt.legend(frameon=False)
plt.savefig("energy_plot.png")


#To get structure dependent pseudo energies we use crambin
crambin = io.LoadPDB("1crn",remote=True).Select("peptide=True")

#Let's get the pseudo energy for the residue num 42 given the full
#strucure as environment
print "local e of res num 42: ", pot.GetEnergy(crambin.residues[41],crambin)

#Let's print out all single residue scores
print "total e of crambin: ", pot.GetTotalEnergy(crambin,crambin)

#We might want all single residue scores...
print "single res energies crambin: ", pot.GetEnergies(crambin,crambin)
