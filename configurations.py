# -*- coding: utf-8 -*-
"""
Created on Mon Jun  1 14:24:30 2015

@author: mathias, aljandra
"""

"""All code for encoding configurations and the operations that can be carried out on them."""

import numpy as np

class configuration_new(object):

    def __init__(self,S,nR,nC):
        validateInput(S,nR,nC)
        self.S = np.matrix(S)
        self.nR = np.array(nR)
        self.nC = np.array(nC)
        self.rows = len(nR)
        self.columns = len(nC)


    def __hash__(self):
        ## x ^ y returns bitwise xor of integers x and y
        #nHash = hash(tuple(i for i in np.sort(self.n)))
#        if len(self.S.shape) < 2:
#        try:
        rowHash    = reduce(lambda x,y: x ^ y, map(hash,(count0123([self.S[i,j] for j in xrange(self.S.shape[1])]) + (self.nR[i],) for i in xrange(self.S.shape[0]))))
        columnHash = reduce(lambda x,y: x ^ y, map(hash,(count0123([self.S[i,j] for i in xrange(self.S.shape[0])]) + (self.nC[j],) for j in xrange(self.S.shape[1]))))
        return rowHash ^ columnHash
#        except:
#            print "WTF!"
#            print self
#            return hash(sum(self.S)) ^ hash(sum(self.n))

    def __str__(self):
#        return "S =\n %s\nn = \n %s"%(str(self.S),str(self.n))
        lines = []
        rowCounts = rowToString(self.nC,sepChar = " ", f = lambda x: str(x))
        rowCountsSep = "="*len(rowCounts)
        lines += [rowCounts, rowCountsSep]
        for i in xrange(self.S.shape[0]):
#            for j in xrange(self.n[i]):
#                lines.append(rowToString([self.S[i,j] for j in xrange(self.S.shape[1])]))
            lines.append(rowToString([self.S[i,j] for j in xrange(self.S.shape[1])]) + " ||  %i"%self.nR[i])
        return "\n".join(lines)

    def __repr__(self):
        return "%s\n(%r)" % (self.__class__, self.__dict__)

    def nullColumns(self):
        #give all indices of columns that are 0
        return [j for j in range(self.L) if max(self.S[:,j]) == 0]

    def segregatingColumns(self):
        return [j for j in range(self.L) if max(self.S[:,j]) != 0]

    def hasRow(self,new_row):
        for row in self.S:
            if np.all(new_row == row):
                return True
        return False

    def whichRow(self,new_row):
        for i,row in enumerate(self.S):
            if np.all(new_row == row):
                return i
        return 0

class configuration(object):
    """
    OLD and depricated version of configuration. Kept in its current state so that
    old code won't be broken. New code relies on new version (which counts rows with
    multiplicity)
    """
    def __init__(self,S,n):
        validateInput(S,n)
        self.S = np.matrix(S)
        self.n = np.array(n)
        self.L = S.shape[1]


    def __hash__(self):
        ## x ^ y returns bitwise xor of integers x and y
        #nHash = hash(tuple(i for i in np.sort(self.n)))
#        if len(self.S.shape) < 2:
#        try:
        extended_rowHash = reduce(lambda x,y: x ^ y , map(hash,(count0123([self.S[i,j] for j in xrange(self.S.shape[1])])+(self.n[i],) for i in xrange(self.S.shape[0]))))
        columnHash = reduce(lambda x,y: x ^ y, map(hash,(count0123([self.S[i,j] for i in xrange(self.S.shape[0])]) for j in xrange(self.S.shape[1]))))
        return extended_rowHash ^ columnHash
#        except:
#            print "WTF!"
#            print self
#            return hash(sum(self.S)) ^ hash(sum(self.n))


    def __str__(self):
#        return "S =\n %s\nn = \n %s"%(str(self.S),str(self.n))
        rowStrings = []
        for i in xrange(self.S.shape[0]):
#            for j in xrange(self.n[i]):
#                rowStrings.append(rowToString([self.S[i,j] for j in xrange(self.S.shape[1])]))
            rowStrings.append(rowToString([self.S[i,j] for j in xrange(self.S.shape[1])]) + " x %i"%self.nR[i])
        return "\n".join(rowStrings)

    def __repr__(self):
        return "%s\n(%r)" % (self.__class__, self.__dict__)

    def nullColumns(self):
        #give all indices of columns that are 0
        return [j for j in range(self.L) if max(self.S[:,j]) == 0]

    def segregatingColumns(self):
        return [j for j in range(self.L) if max(self.S[:,j]) != 0]

    def hasRow(self,new_row):
        for row in self.S:
            if np.all(new_row == row):
                return True
        return False

    def whichRow(self,new_row):
        for i,row in enumerate(self.S):
            if np.all(new_row == row):
                return i
        return 0

#    def sort(self):
#        S_old = np.matrix(self.S)
#        n_old = np.array(self.n)
#        S_sorted_cols = sortColumns(self.S)
#        extendedS = np.c_[S_sorted_cols,self.n]
#        extendedS_sorted = sortRows(extendedS)
#        self.S = extendedS_sorted[:self.S.shape[0],:self.S.shape[1]]
#        self.n = extendedS_sorted[:,-1]
#        if np.any(self.S != S_old) or np.any(self.n != n_old):
#            self.sort()


def validateInput(S,n):
    pass

def validateInput(S,nR,nC):
    if S.shape != nR.shape + nC.shape:
        raise Exception("S.shape != nR.shape + nC.shape")

def count0123(vec):
    # input a collection "vec" of numbers in {0,1,2,3,4}
    #returns a tuple (x0,x1,x2,x3,x4), s.t. xi == count(eintries in x == i)
    return tuple(sum(a==i for a in vec) for i in (0,1,2,3))

def sortRows(A):
    #    newA = np.array(A)
    #    n_columns = A.shape[1]
    #    for i in range(n_columns)[::-1]:
    #        newA = newA[newA[:,i].argsort()[::-1]]
    #    return newA
    rowList = [tuple(row) for row in np.array(A)]
    rowList.sort()
    return np.matrix(rowList)

def sortColumns(A):
    columnList = [tuple(c) for c in np.array(A.T)]
    columnList.sort()
    return np.transpose(np.matrix(columnList))
#    return np.transpose(sortRows(np.transpose(A)))

def rmZeros(phi):
        new_S = phi.S
        new_S = phi.S[:,phi.segregatingColumns()]
        return configuration(S=new_S,n=phi.n)

# def coalesce(phi,i):
#     if phi.n[i] < 2:
#         raise ValueError("cannot coalesce, since n[%i] < 2."%i)
#     else:
#         n_new=phi.n
#         n_new[i]=phi.n[i]-1
#         return configuration(S = phi.S, n=n_new)
def coalesce(phi,i):
    if phi.nR[i] < 2:
        raise ValueError("cannot coalesce, since nR[%i] < 2."%i)
    else:
        n_new=phi.nR
        n_new[i]=phi.nR[i]-1
        return configuration_new(S = phi.S, nR=n_new, nC = phi.nC)


def nucleotideToChar(x):
    if x == 0:
        return "-"
    if x == 1:
        return "+"
    if x == 2:
        return "o"
    if x == 3:
        return "x"

def rowToString(row,width=1,sepChar = " ",f = nucleotideToChar):
    sep = sepChar*width
    return sep+sep.join(map(f,row))+sep

def mutation_new(phi,i,j,X):
    """
    Input:
    phi = (S,nR,nC) -- a configutation
    i,j = row/column index of mutation
    X   = Nucleotide after mutation -- encoded as integer in {0,1,2,3}

    Output:
    phiNew = configuration after applying mutation
    Nr = Combinatorial row-coefficient
    Nc = ocmbinatorial column-coeffiecient
    """

    VERBOSE = False
    # If set to True, the state after each step will be prited.

    if phi.S[i,j] == X:
        return phi,phi.nR[i],phi.nC[j]

    Nr, Nc = 1,1
    nR_new = np.array(tuple(phi.nR)+(1,))
    nR_new[i] = nR_new[i] - 1
    nC_new = np.array(tuple(phi.nC)+(1,))
    nC_new[j] = nC_new[j] - 1

    newRow = np.c_[phi.S[i,:], X]
    S_withExtraCol = np.c_[phi.S,phi.S[:,j]]
    S_new = np.vstack([S_withExtraCol,newRow])

    if VERBOSE:
        print "After step 1"
        print S_new
        print nR_new
        print nC_new
        print ""

    #Check if old rows and columns have multiplicity 0, and if so, delete them
    if nR_new[i] == 0:
        otherRows = [k for k in range(S_new.shape[0]) if k != i]
        S_new = S_new[ otherRows , :]
        nR_new = nR_new[otherRows]

    if VERBOSE:
        print "After step 2.1 (remove rows w multiplicity 0)"
        print S_new
        print nR_new
        print nC_new
        print ""

    if nC_new[j] == 0:
        otherColumns = [l for l in range(S_new.shape[1]) if l != j]
        S_new = S_new[ : , otherColumns]
        nC_new = nC_new[otherColumns]

    if VERBOSE:
        print "After step 2.2 (remove cols w. multiplicity 0)"
        print S_new
        print nR_new
        print nC_new
        print ""

    #Check if new rows/columns already exist, and if so remove them and update mutiplicities
    k = firstRowMatch(S_new,S_new[-1])
    if k < S_new.shape[0] - 1:
        if VERBOSE:
            print "Step 3.1 carried out"
        S_new = S_new[:-1,:]
        nR_new = nR_new[:-1]
        nR_new[k] += 1
        Nr = nR_new[k]

    if VERBOSE:
        print "After step 3.1 (merge identical rows)"
        print S_new
        print nR_new
        print nC_new
        print ""

    l = firstColumnMatch(S_new,S_new[:,-1])
    if l < S_new.shape[1] -1:
        if VERBOSE:
            print "step 3.2 carried out"
            print "l = %i"%l
        S_new = S_new[:,:-1]
        nC_new = nC_new[:-1]
        nC_new[l] += 1
        Nc = nC_new[l]

    if VERBOSE:
        print "After step 3.2 (merge identical columns)"
        print S_new
        print nR_new
        print nC_new
        print ""

    phiNew = configuration_new(S_new,nR_new,nC_new)
    return phiNew, Nr, Nc

def mutation(phi,i,j,k):
    #TODO: Validate input

    #Compute new row
    new_S = np.matrix(phi.S)
    new_S[i,j]=k
    new_row = new_S[i,:]

    #Compare and compute S_new
    if phi.hasRow(new_row):
        #EXISTENT seg. site n[i]=1
        if phi.n[i]<2:
            new_S = np.delete(phi.S,np.s_[i],axis=0) if phi.S.shape[0] > 1 else np.matrix(phi.S)
#            if new_S.shape[0] == 0:
#                print "WTF!!!"
#                print i,j,k
#                print phi
            n_new = np.array(phi.n)
            n_new[phi.whichRow(new_row)] = phi.n[phi.whichRow(new_row)]+1
            n_new = np.delete(n_new,i,0)
            if np.shape(new_S)[1]==1:
                return configuration(S = new_S, n= n_new)
            else:
                 return rmZeros(configuration(S = new_S, n= n_new))
        #EXISTENT seg. site n[i]>=2
        else:
            new_S = phi.S
            n_new = np.array(phi.n)
            n_new[phi.whichRow(new_row)] = phi.n[phi.whichRow(new_row)]+1
            n_new [i] = phi.n[i]-1
            return configuration(S = new_S, n= n_new)
    else:
        #NEW seg. site n[i]=1
        if phi.n[i]<2:
            if np.shape(new_S)[1]==1:
                return configuration(S = new_S, n= phi.n)
            else:
                return rmZeros(configuration(S = new_S, n = phi.n))
        else:
        #EXISTENT seg. site n[i]>=2
            new_S = np.r_[phi.S,new_row]
            n_new = phi.n - np.array([ int(l==i) for l in range(len(phi.n))])
            n_new = np.r_[n_new,1]
            return configuration(S = new_S, n = n_new)

def hasRow(S,new_row):
    for row in S:
        if np.all(new_row == row):
            return True
    return False

def firstRowMatch(S,new_row):
    for i,row in enumerate(S):
        if np.all(new_row == row):
            return i
    return False

def firstColumnMatch(S,new_column):
    return firstRowMatch(np.transpose(S),np.transpose(new_column))

def invisible(phi,i,k):
    #TODO: Validate input

    #n[i]=1
    if phi.n[i]<2:
        extra_col = np.zeros(np.shape(phi.S)[0],dtype=int)
        extra_col[i] = k
        new_S = np.c_[phi.S,extra_col]
        return rmZeros(configuration(S = new_S, n = phi.n))
    else:
        new_S = np.r_[phi.S,phi.S[i,:]]
        extra_col = np.zeros((np.shape(phi.S)[0]+1),dtype=int)
        extra_col[len(extra_col)-1] = k
        new_S = np.c_[new_S,extra_col]
        n_new = np.array(phi.n)
        n_new[i] = phi.n[i]-1
        n_new = np.r_[n_new,1]
        return rmZeros(configuration(S = new_S, n = n_new))

def runTests():
#    print "Hello World"
#
#    A = np.matrix(((1,0,0,1),(1,0,0,0)))
#    b = np.array((1,3),dtype = int)
#
#    phi = configuration(A,b)


#    print "phi.hasRow((1,0,0,1)) = %s"%str(phi.hasRow((1,0,0,1)))
#    print "phi.hasRow((1,0,1,1)) = %s"%str(phi.hasRow((1,0,1,1)))
#    print "phi.whichRow((1,0,0,0)) = %s"%str(phi.whichRow((1,0,0,0)))

    phi = configuration(np.random.randint(0,4,(4,6)),np.array((1,2,1,1)))
    S = phi.S
    for i in range(5):
        print "i=%i"%i
        print "S =\n%s"%str(S)
        S = sortRows(S)
        print "sortRows(S) =\n%s"%str(S)
        S = sortColumns(S)
        print "sortColumns(sortRows(S)) =\n%s\n"%str(S)

#    phi.sort()

#runTests()
