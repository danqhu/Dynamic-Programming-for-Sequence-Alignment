# Seq1 = ['C','C','G','T','T','A','C','G']
# Seq2 = ['C','G','G','T','T','C','G','A']

Seq1 = ['A','A','G','A','T','A']
Seq2 = ['A','A','T','C','T','A','T','A']

# Dynamic Programming Algorithm
#
# Recurrence Function
# D(i,0) = 0; D(0,j) = 0
# For i > 0 and j > 0
# D(i,j) = max[D(i-1,j-1) + sim(S1[i],s2[j]), D(i-1,j), D(i,j-1)]

# Score Matrix and Direction Matrix Calculation
import numpy as np
scoMatrix = np.zeros((len(Seq2)+1,len(Seq1)+1))
dirMatrix = np.zeros((len(Seq2)+1,len(Seq1)+1)) # we define that 0 means ↘︎, 1 means→︎，2 means ↓
dirMatrix[0, 1:] = 1
dirMatrix[1:,0] = 2

# Recurrence function
def recurrence(i, j):
    Dij = {}
    if Seq2[i-1] == Seq1[j-1]:
        Dtl = scoMatrix[i - 1][j - 1] + 1
    else:
        Dtl = scoMatrix[i - 1][j - 1]
    Dl = scoMatrix[i][j-1]
    Dt = scoMatrix[i-1][j]
    tl_l_t_list = [Dtl,Dt,Dl]
    Dij['score']=max(tl_l_t_list)
    Dij['direction']=tl_l_t_list.index(max(tl_l_t_list))

    return Dij

# For D(i,0) = 0; D(0,j) = 0, we know the values in first row and first column of scoMatrix are 0.
# So the calculation is started from the i = 1, j = 1.
for i in range(1, len(Seq2)+1):
    for j in range(1, len(Seq1)+1):
        Dij = recurrence(i,j)
        scoMatrix[i][j] = Dij['score']
        dirMatrix[i][j] = Dij['direction']

# Backtracking
# BTi,BTj = [np.where(scoMatrix == np.max(scoMatrix))[0][0],np.where(scoMatrix == np.max(scoMatrix))[1][0]]
BTi = len(Seq2)
BTj = len(Seq1)
newSeq1 = []
newSeq2 = []

while BTi != 0 and BTj != 0:
    if dirMatrix[BTi,BTj] == 0:
        newSeq1.append(Seq1[BTj-1])
        newSeq2.append(Seq2[BTi-1])
        BTi -= 1
        BTj -= 1
    elif dirMatrix[BTi,BTj] == 1:
        newSeq1.append('-')
        newSeq2.append(Seq2[BTi-1])
        BTi -= 1
        BTj = BTj
    else:
        newSeq1.append(Seq1[BTj-1])
        newSeq2.append('-')
        BTi = BTi
        BTj -= 1
newSeq1.reverse()
newSeq2.reverse()

print(newSeq1)
print(newSeq2)