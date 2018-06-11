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
        extended_rowHash = reduce(lambda x,y: x ^ y , map(hash,(count0123([self.S[i,j] for j in xrange(self.S.shape[1])])+(self.n[i],) for i in xrange(self.S.shape[0]))))
        columnHash = reduce(lambda x,y: x ^ y, map(hash,(count0123([self.S[i,j] for i in xrange(self.S.shape[0])]) for j in xrange(self.S.shape[1]))))
        return extended_rowHash ^ columnHash


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

    def __eq__(self,other):

#        if isinstance(other,node):
#            if self.b != other.b:
#                return False

        if isinstance(other,configuration_new):

            #for bug-checking. Set to false when deploying
            Verbose = False

            #Step 1: verify that row-counts match (up to reordering of rows).
#            extended_rows_self  = [count0123(s) + (self.n[i],) for i,s in enumerate(self.S) ]
#            extended_rows_other = [count0123(s)+(other.n[i],) for i,s in enumerate(other.S)]
            extended_rows_self  = [count0123(tuple(self.S[i,j] for j in range(self.S.shape[1]))) + (self.nR[i],) for i in range(self.S.shape[0]) ]
            extended_rows_other = [count0123(tuple(other.S[i,j] for j in range(other.S.shape[1]))) + (other.nR[i],) for i in range(other.S.shape[0]) ]
            extended_rows_self_sorted  = list(extended_rows_self)
            extended_rows_other_sorted = list(extended_rows_other)
#            if True:
#                print extended_rows_self_sorted
#                print [count0123(s.flat) + (self.n[i],) for i,s in enumerate(self.S) ]
#                print extended_rows_other
            extended_rows_self_sorted.sort()
            extended_rows_other_sorted.sort()


            if not extended_rows_self_sorted == extended_rows_other_sorted:
                if Verbose : print "Rows do not match!"
                return False

            #Step 2: verify that column-counts match up to reordering of columns.
#            extended_columns_self  = [count0123(c) for c in self.S.T ]
#            extended_columns_other = [count0123(c) for c in other.S.T]
            extended_columns_self  = [count0123(tuple(self.S[i,j] for i in xrange(self.S.shape[0]))) + (self.nC[j],) for j in range(self.S.shape[1]) ]
            extended_columns_other = [count0123(tuple(other.S[i,j] for i in xrange(self.S.shape[0]))) + (other.nC[j],) for j in range(self.S.shape[1])]
            extended_columns_self_sorted  = list(extended_columns_self)
            extended_columns_other_sorted = list(extended_columns_other)
            extended_columns_self_sorted.sort()
            extended_columns_other_sorted.sort()

            if not extended_columns_self_sorted == extended_columns_other_sorted:
                if Verbose : print "Columns do not match!"
                return False

            #Step 3: if all row and and column counts match, attempt to
            #        reconstruct. Permutations of rows and columns. Based on
            #        backtracking.

            #work out row and column blocks of self and other.
            # The basic idea is simple. If self and other are equivalent, then
            # any permutation of rows which can turn other into self must
            # satisfy: any element of row_blocks_other[i] is mapped to an
            # element of row_blocks_self[i]
            extended_rows_no_duplicates = list(set(extended_rows_self))
            extended_rows_no_duplicates.sort()
            n_row_blocks = len(extended_rows_no_duplicates)
            row_blocks_self  = [ [] for i in xrange(n_row_blocks)]
            row_blocks_other = [ [] for i in xrange(n_row_blocks)]
            for i in xrange(self.S.shape[0]):
                row_blocks_self[ extended_rows_no_duplicates.index(extended_rows_self[i]) ].append(i)
                row_blocks_other[extended_rows_no_duplicates.index(extended_rows_other[i])].append(i)
            row_blocks_self  = [ tuple(l) for l in row_blocks_self ]
            row_blocks_other = [ tuple(l) for l in row_blocks_other]

            # We now do the same thing for the columns of self and other
            extended_columns_no_duplicates = list(set(extended_columns_self))
            extended_columns_no_duplicates.sort()
            n_column_blocks = len(extended_columns_no_duplicates)
            column_blocks_self  = [ [] for i in xrange(n_column_blocks)]
            column_blocks_other = [ [] for i in xrange(n_column_blocks)]
            for i in xrange(self.S.shape[1]):
                column_blocks_self[ extended_columns_no_duplicates.index(extended_columns_self[i]) ].append(i)
                column_blocks_other[extended_columns_no_duplicates.index(extended_columns_other[i])].append(i)
            column_blocks_self  = [ tuple(l) for l in column_blocks_self ]
            column_blocks_other = [ tuple(l) for l in column_blocks_other]

            S_s = self.S
            S_o = other.S

            if Verbose :
                print row_blocks_self,"\n",row_blocks_other
                print column_blocks_self,"\n",column_blocks_other
                for i in range(len(row_blocks_other)):
                    print i,np.all(S_s[row_blocks_self[i][0]] == S_o[row_blocks_other[i][0]])

            return __reconstruction_possible__(self.S,other.S,row_blocks_self,row_blocks_other,column_blocks_self,column_blocks_other)

#            self.sort()
#            other.sort()
#            return np.all(self.S == other.S) and np.all(self.n == other.n)

        else:
            return NotImplemented

    def __ne__(self,other):
        result = self.__eq__(other)
        if result is NotImplemented:
            return result
        return not result
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
    """
    Count the number of ocurrebnces of 0, 1, 2 and 3 respectively in a vector.

    Args:
       vec -- A vector of integers
    Returns:
       a tuple (x0,x1,x2,x3) such that xi == number of ocurrences of i in vec.
    """
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
    else:
        return "?"

def rowToString(row,width=1,sepChar = " ",f = nucleotideToChar):
    sep = sepChar*width
    return sep+sep.join(map(f,row))+sep

def mutation_new(phi,i,j,X):
    """
    Input:
    phi = (S,nR,nC) -- a configutation
    i,j = row/column index of mutation
    X   = Nucleotide to substitute for S[i,j] -- encoded as integer in {0,1,2,3}

    Output:
    phiNew = configuration after applying mutation
    Nr = Combinatorial row-coefficient (i.e. count of the row produced by subsituting S[i,j] with X, after substitution)
    Nc = ocmbinatorial column-coeffiecient (i.e. count of the column produced by subsituting S[i,j] with X, after substitution)
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

#def sortRows(A):
#    newA = np.array(A)
#    n_columns = A.shape[1]
#    for i in range(n_columns)[::-1]:
#        newA = newA[newA[:,i].argsort()[::-1]]
#    return newA

#def sortColumns(A):
#    return np.transpose(sortRows(np.transpose(A)))

def __reconstruction_possible__(A,
                                B,
                                row_blocks_A = None,
                                row_blocks_B = None,
                                column_blocks_A = None,
                                column_blocks_B = None,
                                sigma_rows_A = (),
                                sigma_rows_B = (),
                                sigma_columns_A = (),
                                sigma_columns_B = ()):
    """
    Test if it is possible to construct a permutation of rows and collumns
    that turn A into B, using a backtracking-algorithm.

    Arguments:
       A,B -- Two matricees satisfying A.shape == B.shape

       row_blocks_A,row_blocks_B -- Two partitions of [0,...,A.shape[0]],
       encoded as a list of tuples.This partiton into blocks is precomputed so
       as to speed up execution of the algorithm, and satisfies that all
       indicees in each block correspond to rows with the same
       nucleotide-counts, i.e. for any i,j and k (within range) the following
       holds:
       count0123(A[row_blocks_A[i][j],:]) == count0123(A[row_blocks_A[i][k],:])
       and
       count0123(B[row_blocks_A[i][j],:]) == count0123(B[row_blocks_B[i][k],:])

       column_blocks_A,column_blocks_B -- Same as row_blocks_A, row_blocks_B,
       only this time a partition of the column-indicees.

       sigma_rows_A,sigma_rows_B,sigma_columns_A,sigma_columns_B -- Rows and
       Columns that are already known to match. The arguments passed to this
       method should always satisfy:
       A[sigma_rows_A,:][:,sigma_columns_A] == B[sigma_rows_B,:][:,sigma_columns_B].

    Returns:
        True if a reconstruction is possible.
    """

    Verbose = False #Produce debugging output if set to true.

    if Verbose:
#        print A[sigma_rows_A,:][:,sigma_columns_A]
#        print B[sigma_columns_B,:][:,sigma_columns_B]
        truthSoFar = A[sigma_rows_A,:][:,sigma_columns_A] == B[sigma_rows_B,:][:,sigma_columns_B]
        if not np.all(truthSoFar):
            print A[sigma_rows_A,:][:,sigma_columns_A] == B[sigma_rows_B,:][:,sigma_columns_B]
            print sigma_rows_A,sigma_columns_A
            print sigma_rows_B,sigma_columns_B


    rowsFilled = len(sigma_rows_A) == A.shape[0]
    columnsFilled = len(sigma_columns_A) == A.shape[1]
    if rowsFilled and columnsFilled:
        assert np.all(A[sigma_rows_A,:][:,sigma_columns_A] == B[sigma_rows_B,:][:,sigma_columns_B])
        return True

    if not rowsFilled:
        # pick blocks of minimal size in row_blocks_A and column_blocks_A
        # this guarantees us that we have to make as few guesses as possible.
        rows_working_block_index= 0
        for i,block in enumerate(row_blocks_A):
            if len(row_blocks_A[rows_working_block_index]) > len(block):
                rows_working_block_index= i

        A_rows_working_block = row_blocks_A[rows_working_block_index]
        B_rows_working_block = row_blocks_B[rows_working_block_index]

        row_blocks_A.remove(A_rows_working_block)
        row_blocks_B.remove(B_rows_working_block)
        if len(A_rows_working_block) > 1:
            row_blocks_A.append(A_rows_working_block[1:])
            appendB_rows_later = True
        else:
            appendB_rows_later = False

        i_A = A_rows_working_block[0]
        new_sigma_rows_A = sigma_rows_A + (i_A,)
    else:
        new_sigma_rows_A = sigma_rows_A

    if not columnsFilled:
        columns_working_block_index = 0
        for i,block in enumerate(column_blocks_A):
            if len(column_blocks_A[columns_working_block_index]) > len(block):
                columns_working_block_index= i

        A_columns_working_block = column_blocks_A[columns_working_block_index]
        B_columns_working_block = column_blocks_B[columns_working_block_index]

        column_blocks_A.remove(A_columns_working_block)
        column_blocks_B.remove(B_columns_working_block)
        if len(A_columns_working_block) > 1:
            column_blocks_A.append(A_columns_working_block[1:])
            appendB_columns_later = True
        else:
            appendB_columns_later = False

        j_A = A_columns_working_block[0]
        new_sigma_columns_A = sigma_columns_A + (j_A,)
    else:
        new_sigma_columns_A = sigma_columns_A


    B_row_column_pairs = []
    if (not columnsFilled) and (not rowsFilled):
        for i in B_rows_working_block:
            possibleRow = sigma_rows_B + (i,)
            for j in B_columns_working_block:
                possibleColumn = sigma_columns_B + (j,)
                if np.all(A[new_sigma_rows_A,:][:,new_sigma_columns_A] == B[possibleRow,:][:,possibleColumn]):
                   B_row_column_pairs.append( (possibleRow , possibleColumn) )


    if columnsFilled and (not rowsFilled):
#        range_i_B = B_rows_working_block
#        test = lambda i: np.all( A[new_sigma_rows_A,:][:,sigma_columns_A] == B[sigma_rows_B + (i,),:][:,sigma_columns_B])
#        matching_i_B = filter(test, B_rows_working_block)
#        B_row_column_pairs = [ (sigma_rows_B + (i,) , sigma_columns_B ) for i in matching_i_B ]
        for i in B_rows_working_block:
            possibleRow = sigma_rows_B + (i,)
            if np.all(A[new_sigma_rows_A,:][:,sigma_columns_A] == B[possibleRow,:][:,sigma_columns_B]):
                B_row_column_pairs.append( (possibleRow , sigma_columns_B) )

    if (not columnsFilled) and rowsFilled:
#       test = lambda j : np.all( A[sigma_rows_A,:][:,new_sigma_columns_A] == B[sigma_rows_B,:][:,sigma_columns_B + (j,) ] )
#       matching_j_B = filter( test , B_columns_working_block )
#       B_row_column_pairs = [ (sigma_rows_B , sigma_columns_B + (j,) ) for j in matching_j_B ]
#       new_B_rows = filter(lambda i: np.all(A[new_sigma_rows_A:][:,sigma_columns_A] == B[sigma_rows_B + (i,),:][:,sigma_columns_B]),B_rows_working_block
        for j in B_columns_working_block:
            possibleColumn = sigma_columns_B + (j,)
            if np.all( A[sigma_rows_A,:][:,new_sigma_columns_A] == B[sigma_rows_B,:][:,possibleColumn] ):
                B_row_column_pairs.append( (sigma_rows_B , possibleColumn) )

#    if len(matching_ij_B) == 0:
#        return False
#    else:
#        for i,j in matching_ij_B:
#            possible_sigma_rows_B = sigma_rows_B + (i,)
#            possible_sigma_columns_B = sigma_columns_B + (j,)
#            if np.all(A[sigma_rows_A,:][:,sigma_columns_A] == B[possible_sigma_rows_B,:][:,possible_sigma_rows_B]):
#
#                if appendB_rows_later:
#                    possible_row_blocks_B = row_blocks_B + [filter(lambda x: x != i,B_rows_working_block)]
#                else:
#                    possible_row_blocks_B = list(row_blocks_B)
#
#                if appendB_columns_later:
#                    possible_column_blocks_B = column_blocks_B + [filter(lambda x: x != i,B_columns_working_block)]
#                else:
#                    possible_column_blocks_B = list(column_blocks_B)
#
#                if reconstruction_possible(
#                    A = A, B = B,
#                    row_blocks_A = list(row_blocks_A), # list(l) passes a COPPY of l instead of l itself
#                    row_blocks_B = possible_row_blocks_B,
#                    column_blocks_A = list(column_blocks_A),
#                    column_blocks_B = possible_column_blocks_B,
#                    sigma_rows_A = sigma_rows_A,
#                    sigma_columns_A = sigma_columns_A,
#                    sigma_rows_B = possible_sigma_rows_B,
#                    sigma_columns_B = possible_sigma_columns_B
#                    ):
#                    return True

    if len(B_row_column_pairs) == 0:
        #The perumtation being explored could not be expanded. we must return false
#        if Verbose:
#            print "No suitable Columns/Rows!"
#            print sigma_rows_A,"--->", new_sigma_rows_A,"\t",sigma_rows_B,"-/->",B_rows_working_block,"\nA[%i,:] = %s\nB[%i,:] = %s"%(i_A,str(A[i_A,:]),B_rows_working_block[0],str(B[B_rows_working_block[0],:]))
#            print "A[new_sigma_rows_A,:][:,new_sigma_columns_A]\n  = %s"%str(A[new_sigma_rows_A,:][:,new_sigma_columns_A])
#            print "B[ sigma_rows_B + (B_rows_working_block[0],) ,:][:,sigma_columns_B + (B_columns_working_block[0],) ]\n  = %s"%str( B[ sigma_rows_B + (B_rows_working_block[0],) ,:][:,sigma_columns_B + (B_columns_working_block[0],) ] )
#            print sigma_columns_A,"-->", new_sigma_columns_A,"\t",sigma_columns_B,"-/->",B_columns_working_block
#            print ""
#            print B_rows_working_block, B_columns_working_block
#            print row_blocks_A
#            print column_blocks_A
#            print row_blocks_B
#            print column_blocks_B
#            print A == B[::-1,:]
        return False
    for new_sigma_rows_B, new_sigma_columns_B in B_row_column_pairs:

        #Build new row- and column-blocks
        if not rowsFilled:
            if appendB_rows_later:
                i = new_sigma_rows_B[-1]
                possible_row_blocks_B = row_blocks_B + [filter(lambda x: x != i, B_rows_working_block)]
            else:
                possible_row_blocks_B = list(row_blocks_B)
        else:
            possible_row_blocks_B = list(row_blocks_B)

        if not columnsFilled:
            if appendB_columns_later:
                j = new_sigma_columns_B[-1]
                possible_column_blocks_B = column_blocks_B + [filter(lambda x: x != j, B_columns_working_block)]
            else:
                possible_column_blocks_B = list(column_blocks_B)
        else:
            possible_column_blocks_B = list(column_blocks_B)
#        pass
#        print new_sigma_rows_A,new_sigma_columns_A
#        print new_sigma_rows_B,new_sigma_columns_B
#        print A[new_sigma_rows_A,:][:,new_sigma_columns_A]
#        print B[new_sigma_rows_B,:][:,new_sigma_columns_B]
#        print np.all(A[new_sigma_rows_A,:][:,new_sigma_columns_A] == B[new_sigma_rows_B,:][:,new_sigma_columns_B])
#        print row_blocks_A
#        print column_blocks_A
#        print row_blocks_B
#        print column_blocks_B
        if __reconstruction_possible__(
                A = A,
                B = B,
                row_blocks_A = list(row_blocks_A), # list(l) passes a COPPY of l instead of l itself
                row_blocks_B = list(possible_row_blocks_B),
                column_blocks_A = list(column_blocks_A),
                column_blocks_B = list(possible_column_blocks_B),
                sigma_rows_A = new_sigma_rows_A,
                sigma_columns_A = new_sigma_columns_A,
                sigma_rows_B = new_sigma_rows_B,
                sigma_columns_B = new_sigma_columns_B
                ):
#                if Verbose: print "Reconstruction Possible!!!"
                return True
    # When all ways of picking permutations have been exhausted, and none have
    # worked, we must return false.
    return False


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
