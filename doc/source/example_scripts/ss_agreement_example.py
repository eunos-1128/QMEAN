from qmean import SSAgreement

#get an SSAgreement object
ss = SSAgreement()

#Let's assume we have a helix as predicted secondary structure.
# => psipred state of H
#We now evaluate what the score is when we find either a
#helix or a sheet in the actual structure with an increasing
#confidence value of the psipred prediction.

print("PSIPRED prediction is H!")

for i in range(10):
    s_one = ss.LogOddScore('H','H',i)
    s_two = ss.LogOddScore('E','H',i)

    print("Actual PSIPRED confidence:",i)
    print("Score given a DSSP state of H:", s_one)
    print("Score given a DSSP state of E:", s_two)
