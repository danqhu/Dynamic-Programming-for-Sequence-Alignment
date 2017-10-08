# -*- coding: utf-8 -*--
# Dynamic Programming Algorithm Global Alignment
#
# Recurrence Function
# D(i,0) = i*Sdel; D(0,j) = j*Sins
# For i > 0 and j > 0
# D(i,j) = max[D(i-1,j-1) + sim(S1[i],s2[j]), D(i-1,j)-Sdel, D(i,j-1)-Sins]
class dynamic_programming_global_alignment():
    def __init__(self, Seq1, Seq2, Ssame, Sdiff, Sdel, Sins):
        self.Seq1 = Seq1
        self.Seq2 = Seq2
        self.Ssame = Ssame
        self.Sdiff = Sdiff
        self.Sdel = Sdel
        self.Sins = Sins

    # Score Matrix and Direction Matrix Calculation
    import numpy as np
    # Initialize the Score Matrix and Direction Matrix
    def __init_score_direction_matrics(self):
        import numpy as np
        self.__scoMatrix = np.zeros((len(Seq2)+1,len(Seq1)+1))
        # self.__scoMatrix[0] = [self.Sdel*i for i in range(len(Seq1)+1)]
        # self.__scoMatrix[:,0] = [self.Sins*j for j in range(len(Seq2)+1)]

        self.__dirMatrix = np.zeros((len(Seq2)+1,len(Seq1)+1)) # we define that 0 means ↘︎, 1 means→︎，2 means ↓
        # self.__dirMatrix[0, 1:] = 1
        # self.__dirMatrix[1:,0] = 2


    # Recurrence function
    def __recurrence(self,i, j):
        Dij = {}
        if j == 0:
            Dij['score'] = self.Sdel*i
            Dij['direction'] = 1
        elif i == 0:
            Dij['score'] = self.Sins * j
            Dij['direction'] = 2
        else:
            if Seq2[i-1] == Seq1[j-1]:
                Dtl = self.__scoMatrix[i - 1][j - 1] + self.Ssame
            else:
                Dtl = self.__scoMatrix[i - 1][j - 1] + self.Sdiff
            Dl = self.__scoMatrix[i][j-1]+self.Sdel
            Dt = self.__scoMatrix[i-1][j]+self.Sins
            tl_l_t_list = [Dtl,Dt,Dl]
            Dij['score']=max(tl_l_t_list)
            Dij['direction']=tl_l_t_list.index(max(tl_l_t_list))

        return Dij

    # For D(i,0) = 0; D(0,j) = 0, we know the values in first row and first column of scoMatrix are 0.
    # So the calculation is started from the i = 1, j = 1.
    def __calculate_score_direction_matrics(self):
        for i in range(0, len(Seq2)+1):
            for j in range(0, len(Seq1)+1):
                Dij = self.__recurrence(i,j)
                self.__scoMatrix[i][j] = Dij['score']
                self.__dirMatrix[i][j] = Dij['direction']
        self.scoMatrix = self.__scoMatrix
        self.dirMatrix = self.__dirMatrix

    # Backtracking
    # BTi,BTj = [np.where(scoMatrix == np.max(scoMatrix))[0][0],np.where(scoMatrix == np.max(scoMatrix))[1][0]]
    def __backtracking(self):

        BTi = len(Seq2)
        BTj = len(Seq1)
        self.newSeq1 = []
        self.newSeq2 = []

        while BTi != 0 and BTj != 0:
            if self.__dirMatrix[BTi,BTj] == 0:
                self.newSeq1.append(Seq1[BTj-1])
                self.newSeq2.append(Seq2[BTi-1])
                BTi -= 1
                BTj -= 1
            elif self.__dirMatrix[BTi,BTj] == 1:
                self.newSeq1.append('-')
                self.newSeq2.append(Seq2[BTi-1])
                BTi -= 1
                BTj = BTj
            else:
                self.newSeq1.append(Seq1[BTj-1])
                self.newSeq2.append('-')
                BTi = BTi
                BTj -= 1
        self.newSeq1.reverse()
        self.newSeq2.reverse()

    def align_pairwise_sequences(self):
        self.__init_score_direction_matrics()
        self.__calculate_score_direction_matrics()
        self.__backtracking()
        return self.newSeq1,self.newSeq2





class dynamic_programming_local_alignment():
    def __init__(self, Seq1, Seq2, Ssame, Sdiff, Sdel, Sins):
        self.Seq1 = Seq1
        self.Seq2 = Seq2
        self.Ssame = Ssame
        self.Sdiff = Sdiff
        self.Sdel = Sdel
        self.Sins = Sins
        # dynamic_programming_global_alignment.__init__(self, Seq1, Seq2, Ssame, Sdiff, Sdel, Sins)

    def __init_score_direction_matrics(self):
        import numpy as np
        self.__scoMatrix = np.zeros((len(Seq2) + 1, len(Seq1) + 1))
        # self.__scoMatrix[0] = [self.Sdel*i for i in range(len(Seq1)+1)]
        # self.__scoMatrix[:,0] = [self.Sins*j for j in range(len(Seq2)+1)]

        self.__dirMatrix = np.zeros((len(Seq2) + 1, len(Seq1) + 1))  # we define that 0 means ↘︎, 1 means→︎，2 means ↓
        # self.__dirMatrix[0, 1:] = 1
        # self.__dirMatrix[1:,0] = 2


    def __recurrence(self,i, j):
        Dij = {}
        if j == 0:
            Dij['score'] = max(0, self.Sdel * i)
            Dij['direction'] = 1
        elif i == 0:
            Dij['score'] = max(0, self.Sins * j)
            Dij['direction'] = 2
        else:
            if Seq2[i - 1] == Seq1[j - 1]:
                Dtl = self.__scoMatrix[i - 1][j - 1] + self.Ssame
            else:
                Dtl = self.__scoMatrix[i - 1][j - 1] + self.Sdiff
            Dl = self.__scoMatrix[i][j - 1] + self.Sdel
            Dt = self.__scoMatrix[i - 1][j] + self.Sins
            tl_l_t_list = [Dtl, Dt, Dl, 0]
            Dij['score'] = max(tl_l_t_list)
            Dij['direction'] = 0 if tl_l_t_list.index(max(tl_l_t_list)) == 3 else tl_l_t_list.index(max(tl_l_t_list))

        return Dij

    def __calculate_score_direction_matrics(self):
        for i in range(0, len(Seq2) + 1):
            for j in range(0, len(Seq1) + 1):
                Dij = self.__recurrence(i, j)
                self.__scoMatrix[i][j] = Dij['score']
                self.__dirMatrix[i][j] = Dij['direction']
        self.scoMatrix = self.__scoMatrix
        self.dirMatrix = self.__dirMatrix

    def __backtracking(self):
        import numpy as np
        BTi,BTj = [np.where(self.__scoMatrix == np.max(self.__scoMatrix))[0][-1],np.where(self.__scoMatrix == np.max(self.__scoMatrix))[1][-1]]
        # BTi = len(Seq2)
        # BTj = len(Seq1)
        self.newSeq1 = []
        self.newSeq2 = []

        while (BTi != 0 and BTj != 0) and (self.__scoMatrix[BTi,BTj] != 0):
            if self.__dirMatrix[BTi, BTj] == 0:
                self.newSeq1.append(Seq1[BTj - 1])
                self.newSeq2.append(Seq2[BTi - 1])
                BTi -= 1
                BTj -= 1
            elif self.__dirMatrix[BTi, BTj] == 1:
                self.newSeq1.append('-')
                self.newSeq2.append(Seq2[BTi - 1])
                BTi -= 1
                BTj = BTj
            else:
                self.newSeq1.append(Seq1[BTj - 1])
                self.newSeq2.append('-')
                BTi = BTi
                BTj -= 1
        self.newSeq1.reverse()
        self.newSeq2.reverse()


    def align_pairwise_sequences(self):
        self.__init_score_direction_matrics()
        self.__calculate_score_direction_matrics()
        self.__backtracking()
        return self.newSeq1,self.newSeq2







# # Test for class dynamic_programming_global_alignment
#
# # For parameters Ssame = 1, Sdiff = 0, Sdel = 0, Sins = 0
# Seq1 = ['A','A','G','A','T','A']
# Seq2 = ['A','A','T','C','T','A','T','A']
#
# test = dynamic_programming_global_alignment(Seq1,Seq2,1,-1,-1,-1)
# print(test.align_pairwise_sequences())
# print(test.newSeq1)
# print(test.newSeq2)
# print(test.scoMatrix)
# print(test.dirMatrix)
#
# # For parameter Ssame = 1, Sdiff = -1, Sdel = -1, Sins = -1
# Seq1 = ['G','C','G','T','A','G']
# Seq2 = ['A','C','G','A','T','A','C','C']
#
# test = dynamic_programming_global_alignment(Seq1,Seq2,1,-1,-1,-1)
# print(test.align_pairwise_sequences())
# print(test.newSeq1)
# print(test.newSeq2)
# print(test.scoMatrix)
# print(test.dirMatrix)

# Test for class dynamic_programming_local_alignment
# For parameter Ssame = 1, Sdiff = -1, Sdel = -1, Sins = -1

Seq1 = ['G','C','G','T','A','G']
Seq2 = ['A','C','G','A','T','A','C','C']

testforlocal = dynamic_programming_local_alignment(Seq1,Seq2,1,-1,-1,-1)
print(testforlocal.align_pairwise_sequences())


