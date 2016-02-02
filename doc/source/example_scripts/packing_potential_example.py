#Example code to demonstrate the usage of setting up a packing potential
#and get some basic info out of it. You need matplotlib to run this script!

#load required modules
import os
from qmean import PackingStatistic, PackingPotential, \
                  ChemType
import matplotlib.pyplot as plt


#define some training targets used later on
training_targets = ["7ODC","2QUD","3T47","4C1A","3ZSJ","4I4T","3H6R","3BSO",
                    "3M0Z","2GPI","3MHP","3V19","2IXS","4G6T","3ZG9","3WSG",
                    "1JH6","2XZE","2OLR","1TIF","3C4S","2CVI","4EET","1LM5",
                    "1MJ5","2FJR","2D3G","3ZNV","3WA2","3WU2","4M7R","2PR7",
                    "3FTJ","1KD8","1HBN","4TRK","2GB4","3HNY","4R7Q","1EAQ",
                    "2O1Q","4DX5","1XAK","5CSM","2XWP","2UWA","3SQZ","4MT8",
                    "3HTU"]


#Create a packing statistic considering other atoms within 5 Angstrom
#and counts up to 32
stat = PackingStatistic(5.0,32)


#Let's fill the stat with some data!
for t in training_targets:
    prot_path = os.path.join("example_data",t+".pdb")
    prot = io.LoadPDB(prot_path).Select("peptide=true")
    stat.Extract(prot,prot)


#Create a potential object from the extracted statistic
#We use the uniform reference state for that...
pot = PackingPotential.Create(stat,reference_state="uniform")

#Save the statistic and potential objects
stat_out_path = os.path.join("example_out","packing_statistic.dat")
pot_out_path = os.path.join("example_out","packing_potential.dat")
stat.Save(stat_out_path)
pot.Save(pot_out_path)


#We're not done yet! Let's play around a bit

#Let's investigate the packing behaviour of two atoms in Lysine. 
#The hypothesis is, that the CA tends to have ist energetic minimum
#at more atoms within the cutoff than the NZ.

ca_pseudo_energies = list()
nz_pseudo_energies = list()

for i in range(33):
    count = i
    e = pot.GetEnergy(ChemType.C_LYS_A, count)
    ca_pseudo_energies.append(e)
    e = pot.GetEnergy(ChemType.N_LYS_Z, count)
    nz_pseudo_energies.append(e)

#let's produce an awesome plot and save it down
x = list()
y1 = list()
y2 = list()

for i in range(33):
    x.append(i)
    x.append(i+1)
    y1.append(ca_pseudo_energies[i])
    y1.append(ca_pseudo_energies[i])
    y2.append(nz_pseudo_energies[i])
    y2.append(nz_pseudo_energies[i])

plt.plot(x,y1,label="LYS-CA",color="#004586")
plt.plot(x,y2,label="LYS-NZ",color="#ff420e")
plt.plot([0,33],[0,0],'k')
plt.xlabel("other atoms within cutoff")
plt.ylabel("pseudo energy")
plt.legend(frameon=False)
out_path = os.path.join("example_out","packing_energy_plot.png")
plt.savefig(out_path)
