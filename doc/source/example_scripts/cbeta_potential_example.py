#Example code to demonstrate the usage of setting up a cbeta potential
#and get some basic info out of it in comparison to an interaction potential. 
#You need matplotlib to run this script!

#load required modules
from qmean import CBetaStatistic, CBetaPotential
from qmean import InteractionStatistic, InteractionPotential
from qmean import ChemType
import os
from ost import conop
import matplotlib.pyplot as plt


#define some training targets used later on
training_targets = ["7ODC","2QUD","3T47","4C1A","3ZSJ","4I4T","3H6R","3BSO",
                    "3M0Z","2GPI","3MHP","3V19","2IXS","4G6T","3ZG9","3WSG",
                    "1JH6","2XZE","2OLR","1TIF","3C4S","2CVI","4EET","1LM5",
                    "1MJ5","2FJR","2D3G","3ZNV","3WA2","3WU2","4M7R","2PR7",
                    "3FTJ","1KD8","1HBN","4TRK","2GB4","3HNY","4R7Q","1EAQ",
                    "2O1Q","4DX5","1XAK","5CSM","2XWP","2UWA","3SQZ","4MT8",
                    "3HTU"]


#Create a cbeta statistic, considering pairwise
#distances from 0.0 to 10 A with a binsize of 1 A.
#Pairwise interactions from residues separated by less
#than 4 residues are neglected.
cbeta_stat = CBetaStatistic(0.0,10.0,10,4)

#Create an interaction statistic, considering pairwise
#distances from 0.0 to 10 A with a binsize of 1 A.
#Pairwise interactions from residues separated by less
#than 4 residues are neglected.
interaction_stat = InteractionStatistic(0.0,10.0,10,4)


#Let's fill the stats with some data!
for t in training_targets:
    prot_path = os.path.join("example_data",t+".pdb")
    prot = io.LoadPDB(prot_path).Select("peptide=true")
    cbeta_stat.Extract(prot,prot)
    interaction_stat.Extract(prot,prot)


#Create potential objects from the extracted statistics
cbeta_pot = CBetaPotential.Create(cbeta_stat)
interaction_pot = InteractionPotential.Create(interaction_stat)


#We're not done yet! Let's play around a bit
#Let's investigate the interaction between the CB atoms of Alanine
#as they are stored in the cbeta potential, compared with exactly the
#same atoms in the interaction potential.

pseudo_cbeta_energies = list()
pseudo_interaction_energies = list()

for i in range(10):
    distance = i

    #the cbeta energy can be extracted with the amino acid describtion from ost
    e = cbeta_pot.GetEnergy(conop.ALA, conop.ALA, distance)
    pseudo_cbeta_energies.append(e)

    #The interaction energy has to be extracted with the ChemType
    e = interaction_pot.GetEnergy(ChemType.C_ALA_B, ChemType.C_ALA_B, distance)
    pseudo_interaction_energies.append(e)


#let's produce an awesome plot and save it down
x = list()
y_cbeta = list()
y_interaction = list()

for i in range(10):
    x.append(i)
    x.append(i+1)
    y_cbeta.append(pseudo_cbeta_energies[i])
    y_cbeta.append(pseudo_cbeta_energies[i])
    y_interaction.append(pseudo_interaction_energies[i])
    y_interaction.append(pseudo_interaction_energies[i])

plt.plot(x,y_cbeta,label="ALA-CB  ALA-CB cbeta pot",color='#004586')
plt.plot(x,y_interaction,label="ALA-CB  ALA-CB interaction pot",color="#ff420e")
plt.plot([0,10],[0,0],'k')
plt.xlabel("distance")
plt.ylabel("pseudo energy")
plt.legend(frameon=False)
fig_path = os.path.join("example_out","cbeta_energy_plot.png")
plt.savefig(fig_path)
