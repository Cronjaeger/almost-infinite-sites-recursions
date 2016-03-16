# -*- coding: utf-8 -*-
"""
Created on Mon Jun  1 17:01:15 2015

@author: mathias, alejandra
"""

from configurations import configuration, configuration_new, sortRows, sortColumns
import numpy as np

class node(configuration_new):

#    def __cmp__(self,other):
#        #returnrn 1 if self > other, 0 if self==other, -1 if self < other
##        sum_cmp = cmp(sum(self.n),sum(other.n))
##        if sum_cmp != 0:if np.shape(new_S)[1]==1:
##            return sum_cmp
##        else:
##            if len(s)
#        return cmp(str(self.n),str(other.n)) or cmp(str(self.S),str(other.S))

#    def __init__(self,S,n,b):
#        validateInput(S,n,b)
#        self.S = S
#        self.n = n
#        self.L = S.shape[1]
#        self.b = b

    def __eq__(self,other):

#        if isinstance(other,node):
#            if self.b != other.b:
#                return False

        if isinstance(other,node) or isinstance(other,configuration_new):

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

            return reconstruction_possible(self.S,other.S,row_blocks_self,row_blocks_other,column_blocks_self,column_blocks_other)

#            self.sort()
#            other.sort()
#            return np.all(self.S == other.S) and np.all(self.n == other.n)

        else:
            return NotImplemented


###
# OLD VERSION BELOW
###
#         if isinstance(other,node) or isinstance(other,configuration):
#
#             #for bug-checking. Set to false when deploying
#             Verbose = False
#
#             #Step 1: verify that row-counts match (up to reordering of rows).
# #            extended_rows_self  = [count0123(s) + (self.n[i],) for i,s in enumerate(self.S) ]
# #            extended_rows_other = [count0123(s)+(other.n[i],) for i,s in enumerate(other.S)]
#             extended_rows_self  = [count0123(tuple(self.S[i,j] for j in range(self.S.shape[1]))) + (self.n[i],) for i in range(self.S.shape[0]) ]
#             extended_rows_other = [count0123(tuple(other.S[i,j] for j in range(other.S.shape[1]))) + (other.n[i],) for i in range(other.S.shape[0]) ]
#             extended_rows_self_sorted  = list(extended_rows_self)
#             extended_rows_other_sorted = list(extended_rows_other)
# #            if True:
# #                print extended_rows_self_sorted
# #                print [count0123(s.flat) + (self.n[i],) for i,s in enumerate(self.S) ]
# #                print extended_rows_other
#             extended_rows_self_sorted.sort()
#             extended_rows_other_sorted.sort()
#
#
#             if not extended_rows_self_sorted == extended_rows_other_sorted:
#                 if Verbose : print "Rows do not match!"
#                 return False
#
#             #Step 2: verify that column-counts match up to reordering of columns.
# #            extended_columns_self  = [count0123(c) for c in self.S.T ]
# #            extended_columns_other = [count0123(c) for c in other.S.T]
#             extended_columns_self  = [count0123(tuple(self.S[i,j] for i in xrange(self.S.shape[0]))) for j in range(self.S.shape[1]) ]
#             extended_columns_other = [count0123(tuple(other.S[i,j] for i in xrange(self.S.shape[0]))) for j in range(self.S.shape[1])]
#             extended_columns_self_sorted  = list(extended_columns_self)
#             extended_columns_other_sorted = list(extended_columns_other)
#             extended_columns_self_sorted.sort()
#             extended_columns_other_sorted.sort()
#
#             if not extended_columns_self_sorted == extended_columns_other_sorted:
#                 if Verbose : print "Columns do not match!"
#                 return False
#
#             #Step 3: if all row and and column counts match, attempt to
#             #        reconstruct. Permutations of rows and columns. Based on
#             #        backtracking.
#
#             #work out row and column blocks of self and other.
#             # The basic idea is simple. If self and other are equivalent, then
#             # any permutation of rows which can turn other into self must
#             # satisfy: any element of row_blocks_other[i] is mapped to an
#             # element of row_blocks_self[i]
#             extended_rows_no_duplicates = list(set(extended_rows_self))
#             extended_rows_no_duplicates.sort()
#             n_row_blocks = len(extended_rows_no_duplicates)
#             row_blocks_self  = [ [] for i in xrange(n_row_blocks)]
#             row_blocks_other = [ [] for i in xrange(n_row_blocks)]
#             for i in xrange(self.S.shape[0]):
#                 row_blocks_self[ extended_rows_no_duplicates.index(extended_rows_self[i]) ].append(i)
#                 row_blocks_other[extended_rows_no_duplicates.index(extended_rows_other[i])].append(i)
#             row_blocks_self  = [ tuple(l) for l in row_blocks_self ]
#             row_blocks_other = [ tuple(l) for l in row_blocks_other]
#
#             # We now do the same thing for the columns of self and other
#             extended_columns_no_duplicates = list(set(extended_columns_self))
#             extended_columns_no_duplicates.sort()
#             n_column_blocks = len(extended_columns_no_duplicates)
#             column_blocks_self  = [ [] for i in xrange(n_column_blocks)]
#             column_blocks_other = [ [] for i in xrange(n_column_blocks)]
#             for i in xrange(self.S.shape[1]):
#                 column_blocks_self[ extended_columns_no_duplicates.index(extended_columns_self[i]) ].append(i)
#                 column_blocks_other[extended_columns_no_duplicates.index(extended_columns_other[i])].append(i)
#             column_blocks_self  = [ tuple(l) for l in column_blocks_self ]
#             column_blocks_other = [ tuple(l) for l in column_blocks_other]
#
#             S_s = self.S
#             S_o = other.S
#
#             if Verbose :
#                 print row_blocks_self,"\n",row_blocks_other
#                 print column_blocks_self,"\n",column_blocks_other
#                 for i in range(len(row_blocks_other)):
#                     print i,np.all(S_s[row_blocks_self[i][0]] == S_o[row_blocks_other[i][0]])
#
#             return reconstruction_possible(self.S,other.S,row_blocks_self,row_blocks_other,column_blocks_self,column_blocks_other)
#
# #            self.sort()
# #            other.sort()
# #            return np.all(self.S == other.S) and np.all(self.n == other.n)
#
#         else:
#             return NotImplemented

    def __ne__(self,other):
        result = self.__eq__(other)
        if result is NotImplemented:
            return result
        return not result

#    def __hash__(self):
#        ## x ^ y returns bitwise xor of integers x and y
##        nHash = hash(tuple(i for i in np.sort(self.n)))
#        extended_rowHash = reduce(lambda x,y: x ^ y , map(hash,(count0123(self.S[i,:])+(self.n[i],) for i in xrange(self.S.shape[0]))))
#        columnHash = reduce(lambda x,y: x ^ y, map(hash,(count0123(self.S[:,j]) for j in xrange(self.S.shape[1]))))
#        return extended_rowHash ^ columnHash
#
#    def sort(self):
#        S_sorted_cols = sortColumns(self.S)
#        extendedS = np.c_[S_sorted_cols,self.n]
#        extendedS_sorted = sortRows(extendedS)
#        self.S = extendedS_sorted[:self.S.shape[0],:self.S.shape[1]]
#        self.n = extendedS_sorted[:,-1]

#    def __str__(self):
#        return "S =\n%s\nn =\n%s"%(str(self.S),str(self.n))
#
#    def __repr__(self):
#        return "%s\n(%r)" % (self.__class__, self.__dict__)

class configTable(object):
    """
    A hash-table for storing the results of computations.
    """

    def __init__(self):

        self.__node_table__ = dict()
        self.__size__ = 0
        self.__size_expanded__ = 0
#        self.__edge_count__ = 0
#        self.__b_set__ = set()

    def contains(self,phi,b = None):
        """
        Test if phi is a key in __node_table__.
        If b != None, this method will only return true, if phi is in
        __node__table__ AND a probability is associated with the pair (phi,b)
        """
        if phi not in self.__node_table__:
            return False
        if b == None:
            return True
        else:
            return b in self.__node_table__[phi]

    def get_p(self,phi,b):
        return self.__node_table__[phi][b]

#    def get_m(self,phi):
#        return self.__node_table__[phi][1]

    def add(self,phi,p,b):
#        phiNode = node(phi.S, phi.n)
        phiNode = node(phi.S, phi.nR, phi.nC)
        if phi not in self.__node_table__:
#            newNode = node(phi.S, phi.n)
            self.__node_table__[phiNode] = dict()
            self.__size__ += 1
        self.__node_table__[phiNode][b] = p
        self.__size_expanded__ += 1

#    def set_m(self,phi,m):
##        phiAsNode = node(phi.S,phi.n)
##        self.nodes[phiAsNode][1] = m
#        self.__node_table__[phi][1] = m

#    def set_p(self,phi,p):
#        self.__node_table__[phi] = p

    def get_all_configurations(self):
        """
        if X is an instance of configTable, calling X.get_all_configurations()
        will return a dictionary of configurations mapped to probabilities
        and m-values
        """
        return dict(self.__node_table__)

    def get_size(self):
        """
        returns self.size
        """
        return self.__size__,self.__size_expanded__

#def count0123more(vec):
#    # input a collection "vec" of numbers in {0,1,2,3,4}
#    #returns a tuple (x0,x1,x2,x3,x4), s.t. xi == count(eintries in x leq i)
#    return tuple(sum(a<=i for a in vec) for i in (0,1,2,3,4))
#
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

#def sortRows(A):
#    newA = np.array(A)
#    n_columns = A.shape[1]
#    for i in range(n_columns)[::-1]:
#        newA = newA[newA[:,i].argsort()[::-1]]
#    return newA

#def sortColumns(A):
#    return np.transpose(sortRows(np.transpose(A)))
def reconstruction_possible(A,
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
       as to speed up execution of the algorithm, and satisfies the following:
       For any i,j and k,
       count0123(A[row_blocks_A[i][j],:]) == count0123(A[row_blocks_A[i][j],:])
       must be true.

       column_blocks_A,column_blocks_B -- Same as row_blocks_A, row_blocks_B,
       only this time a partition of the column-indicees.

       sigma_rows_A,sigma_rows_B,sigma_columns_A,sigma_columns_B -- Rows and
       Columns that are already known to match. The arguments passed to this
       method should always satisfy:
       A[sigma_rows_A,:][:,sigma_columns_A] == B[sigma_rows_B,:][:,sigma_columns_B].

    Returns:
        True if a reconstruction is possible.
    """

#    cmp_cardinality = lambda x,y: cmp(len(x),len(y))
#    row_blocks_A.sort(cmp = cmp_cardinality)
#    row_blocks_B.sort(cmp = cmp_cardinality)
#    column_blocks_A.sort(cmp = cmp_cardinality)
#    column_blocks_B.sort(cmp = cmp_cardinality)

#    if sigma_rows_A != None: # verify that things work so far
#        if np.any(A[sigma_rows_A,:][:,sigma_columns_A] != B[sigma_rows_B,:][:,sigma_rows_B]):
#                return False

#    if row_blocks_A == None or row_blocks_B == None or column_blocks_A == None or column_blocks_B == None:
#        row_blocks_A, row_blocks_B, column_blocks_A, column_blocks_B = calculateBlocks(A,B)
#    elif row_blocks_A == [] and row_blocks_B == [] and column_blocks_A == [] and column_blocks_B == []:

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
#        A_ij = A[i_A,j_A]
#        matching_ij_B = [(i,j) for i in B_rows_working_block for j in B_columns_working_block if B[i,j] == A_ij]
#        range_ij_B = [(i,j) for i in B_rows_working_block for j in B_columns_working_block]
#        test = lambda t: np.all( A[new_sigma_rows_A,:][:,new_sigma_columns_A] == B[ sigma_rows_B + (t[0],),:][:, sigma_rows_B + (t[1],)])
#        matching_ij_B = filter( test , range_ij_B )
#        B_row_column_pairs = [ (sigma_rows_B + (t[0],) , sigma_columns_B + (t[1],) ) for t in matching_ij_B ]
#        if Verbose:
#            print "B_row_column_pairs",B_row_column_pairs
#            print "[(sigma_rows_A, sigma_columns_B)]",[(sigma_rows_A, sigma_columns_B)]
#            print [A[new_sigma_rows_A,:][:,new_sigma_columns_A] == B[t[0],:][:,t[1]] for t in B_row_column_pairs]
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
        #The perumtation being expored could not be expanded. we must return false
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
        if reconstruction_possible(
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

# def calculateBlocks(A,B):
#     return [],[],[],[] # TODO: implement this

def runTest():
    """
    Runs a variety of tests to verify that the program performs as intended
    (OUTDATED)
    """
    pass
# #    phi1 = configuration(np.random.randint(0,4,(2,3)),np.array((1,2)))
#     phiList = []
#     phiList_mod = []
#     for  i in range(10):
#         phi = configuration(np.random.randint(0,4,(10,10)),np.repeat(range(1,3),5))
#         phi_ReverseRows = configuration(phi.S[::-1,:],phi.n[::-1])
#         phiList.append(phi)
#         phiList_mod.append(phi_ReverseRows)

#    node1,node2 = node(phi.S,phi.n),node(phi.S[[0,1],:],phi.n[[0,1]])
#    print phi
#    print node1
#    print node2
#    from copy import deepcopy
#    node3 = deepcopy(node1)
#
#    if node3.S[0,0] == 0:
#        node3.S[0,0] = 1
#    else:
#        node3.S[0,0] = 0

#    phi2 = deepcopy(phi)

#    print hash(phi),hash(phi2),hash(node1),hash(node2),hash(node3)
#    print phi

#    phi.sort()
#    for i in range(200):
#        phi.sort()
#        phi2.sort()
#        print np.all(phi2.S == phi.S),np.all(phi2.n == phi.n),hash(phi)==hash(phi2)
##    print hash(phi)
#    print phi
#    phi.sort()
#    print hash(phi)
#    print phi

#     b = 0
#     myTable = configTable()
#     for phi in phiList:
#         myTable.add(phi,p=0.4,b=0)
# #        if myTable.contains(phi):
# #            print phi.S
# #            print "\n"
# #        print phi
#
#     print [myTable.contains(phi,b) for phi in myTable.get_all_configurations()]
#     print [myTable.contains(phi,b) for phi in phiList]
#     print [myTable.contains(phiR,b) for phiR in phiList_mod]
#
#     pathPi = node(S = np.matrix("1,1,0,0;1,0,0,0;0,0,1,1"), n = np.array((1,1,2)))
#     pathPi2 = node(S = np.matrix("1,1,0,0;0,0,1,0;0,0,1,1"), n = np.array((1,1,2)))
#     print pathPi
#     print pathPi2
#     print "hash(pathPi)-hash(pathPi2) =",hash(pathPi)-hash(pathPi2)
#     print "(pathPi == pathPi2) ==",pathPi == pathPi2

#    print map(hash,(count0123([pathPi.S[i,j] for j in xrange(3)])+(pathPi.n[i],) for i in xrange(pathPi.S.shape[0])))
#    print map(hash,(count0123([pathPi2.S[i,j] for j in xrange(3)])+(pathPi2.n[i],) for i in xrange(pathPi2.S.shape[0])))
#    print [count0123(tuple(pathPi.S[i,:]))+(pathPi.n[i],) for i in xrange(pathPi.S.shape[0])]

#    print "\n"
#    for phi in  myTable.__node_table__.keys():
#        print phi.S

#    phi1 = node(S = np.matrix("0,0,1;1,1,0"), n = np.array((1,2)))
#    myTable.add(phi1)
#    print phi1.S
#    print myTable.get_all_configurations().keys()[0].S
#    print np.all(myTable.get_all_configurations().keys()[0].S == phi1.S)
##    print myTable.get_all_configurations().keys()[0] == phi1
#    firstKey = myTable.get_all_configurations().keys()[0]
#    print type(phi1)
#    print type(firstKey)
#    print phi1 == firstKey


#    print myTable.contains(phi)
#    print phi
#    print phi == myTable.__node_table__.keys()[0]

#    print hash(myTable.__node_table__.keys()[0]),hash(phi)
#    print myTable.__node_table__.keys()[0] == node(phi.S,phi.n)
#    print myTable.contains(phi)

#    print hash(phi)
#    phi.sort()
#    print myTable.contains(phi)
#    print myTable.get_p(phi)
#    print myTable.get_m(phi)
#    myTable.set_p(phi,0.5)
#    myTable.set_m(phi,3)
#    print myTable.get_p(phi)
#    print myTable.get_m(phi)

#    print myTable.get_p(phi2)

#    d = {node1:0.1 , node3:0.9}
#    print "d[node1] = %f"%d[node1]
#    print "d[node2] = %f"%d[node2]
#    print "d[node3] = %f"%d[node3]

if __name__ == '__main__':
    runTest()
