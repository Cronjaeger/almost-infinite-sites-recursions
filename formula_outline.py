# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 14:40:17 2015

@author: mathias & alejandra
"""

from configurationTable import configTable,node
from configurations import configuration,configuration_new,coalesce,invisible,mutation,mutation_new
import numpy as np
import probabilities as proba
from copy import deepcopy
import time
import sys


phiTable_std = configTable()
#phi_null = node( S = np.matrix((0,)) , nR = np.array((1,)) )
def stdBoundaryCondition(phi):
    if np.all(phi.S == 0) and phi.S.shape == (1,1):
        return 1.0
    else:
        return 0

#for b in xrange(10):
#    phiTable.add(phi_null,p=1.0,b=0)

def prob_External(S,nR,nC, b, theta = 1.0, returnTable = True ,P = (np.ones((4,4)) - np.eye(4))/3):
    myPhiTable = configTable()
    # phi_null = node( S = np.matrix((0,)) , n = np.array((1,)) )
    # renewPhiTable()
    phi = configuration_new(S,nR,nC)
    T = sum(nC)
    p = prob(phi,theta,b,T,P, myPhiTable)
    size = myPhiTable.get_size()
    if returnTable:
        return p,myPhiTable
    else:
        return p

    # myPhiTable = configTable()
    # # phi_null = node( S = np.matrix((0,)) , n = np.array((1,)) )
    # # renewPhiTable()
    # phi = configuration(S,n)
    # p = prob(phi,theta,b,L,P, myPhiTable)
    # size = myPhiTable.get_size()
    # return p,size

def prob(phi,theta,b,T,P,phiTable = phiTable_std,boundaryCondition = stdBoundaryCondition):

    VERBOSE = False
    # if VERBOSE:
    #     print "considering configuration:"
    #     print phi

    if phiTable.contains(phi,b):
        return phiTable.get_p(phi,b)

    # if phi == phi_null:
    #     phiTable.add( phi = phi , p = 1.0 , b = b)
    #     return 1.0
    if sum(phi.nR) == 1:
        p = boundaryCondition(phi)
        phiTable.add(phi = phi, p = p, b = b)
        return p

    # Compute the minimal number of muattions required
    # minMutations = phi.S.shape[1]
    # deviants = set((1,2,3))
    # minMutations = sum([len(deviants.intersection(set([phi.S[i,j] for i in xrange(phi.S.shape[0])]))) for j in xrange(phi.S.shape[1]) ])
    if boundaryCondition == stdBoundaryCondition:
        deviants = set((1,2,3))
        minMutations = sum( [ len(deviants.intersection(set([phi.S[i,j] for i in xrange(phi.S.shape[0])]))) * phi.nC[j] for j in xrange(phi.S.shape[1]) ] )
    else:
        deviants = set((0,1,2,3))
        minMutations = sum( [(len(deviants.intersection(set([phi.S[i,j] for i in xrange(phi.S.shape[0])]))) - 1 ) * phi.nC[j] for j in xrange(phi.S.shape[1]) ] )


    if minMutations > b:
        if VERBOSE:
            print "Configuration:"
            print phi
            print "has (b,minMutations) = (%i,%i)\n"%(b,minMutations)
        phiTable.add( phi = phi , p = 0.0 , b = b)
        return 0.0

    coalescentIndicees = (i for i in range(len(phi.nR)) if phi.nR[i]>1)
    #mutationIndicees = ((i,j,k) for i in xrange(phi.S.shape[0]) for j in xrange(phi.S.shape[1]) for k in xrange(4) if k != phi.S[i,j] )
    mutationIndicees = ( (i,j,Y) for i in xrange(phi.S.shape[0]) for j in xrange(phi.S.shape[1]) for Y in xrange(4) if Y != phi.S[i,j] )

    fac_c,pro_c,fac_i,pro_i,fac_m,pro_m = 0,0,0,0,0,0

    line_1 = 0
    for i in coalescentIndicees:
        phi_coal = deepcopy(phi)
        fac_c = proba.p_coalesce_new(phi_coal,i,theta)
        pro_c = prob( coalesce(phi_coal,i), theta , b , T, P, phiTable)
        line_1 += fac_c * pro_c

    line_2 = 0
    if b > 0:
        for i,j,Y in mutationIndicees:
            phi_mut,Nr,Nc = mutation_new(deepcopy(phi),i,j,Y)
            fac_m = proba.p_mutation_new(phi, theta, P, X = phi.S[i,j], Y = Y, Nr = Nr, Nc = Nc)
            pro_m = prob( phi_mut, theta, b-1, T, P, phiTable)
            line_2 += fac_m * pro_m

    p_stateChangeOccurs = 1.0
    if P[0,0] + P[1,1] + P[2,2] + P[3,3] > 0.0:
        p_noStatechange = 0.0
        n_phi = sum(phi.nR)
        #T_phi = T
        T_phi = sum(phi.nC)
        invariantFactor = float(theta)/((n_phi - 1 + theta) * T_phi * n_phi)
        for i in xrange(phi.S.shape[0]):
            Nr = phi.nR[i]
            for j in xrange(phi.S.shape[1]):
                Nc = phi.nC[j]
                X = phi.S[i,j]
                p_noStatechange += invariantFactor * Nr * Nc * P[X,X]
        p_stateChangeOccurs -= p_noStatechange

#     line_2 = 0
#     if b > minMutations + 1:
#         invisibleIndicees = ((i,k) for k in xrange(1,4) for i in xrange(phi.S.shape[0]))
#         for i,k in invisibleIndicees:
#             phi_inv = deepcopy(phi)
#             fac_i = proba.p_invisible(phi_inv,i,k,theta,T,P)
#             pro_i = prob( invisible(phi_inv,i,k) , theta , b - 1, T, P, phiTable)
#             line_2 += fac_i * pro_i
#
#     line_3 = 0
#     if b > 0:
#         for i,j,k in mutationIndicees:
#             if ((np.all(phi.S==0))==False):
#                 phi_mut = deepcopy(phi)
#                 fac_m = proba.p_mutation(phi_mut,i,j,k,theta,T,P)
#                 pro_m = prob( mutation(phi_mut,i,j,k), theta , b - 1, T, P, phiTable )
#                 line_3 += fac_m * pro_m
#
# #    print ""
# #    print phi,b
# #    print line_1,fac_c,pro_c,len(list(coalescentIndicees))
# #    print line_2,fac_i,pro_i,len(list(invisibleIndicees))
# #    print line_3,fac_m,pro_m,len(list(mutationIndicees))
#     p = line_1 + line_2 + line_3
#     phiTable.add( phi = phi , p = p , b = b)
# #    print "womp"
    p = (line_1 + line_2)/p_stateChangeOccurs
    phiTable.add(phi = phi, p = p, b = b)
    return p

if __name__ == "__main__":

    sys.setrecursionlimit(1000000)

    P = (np.ones((4,4)) - np.eye(4))/3
    theta = np.array((0.14,14,1400),float)
    T = 14000

    #Easy data set
    S_1 = np.matrix((0,))
    n_1 = np.array((5,))
    phi_1 = configuration(S_1,n_1)

    #Typical example
    S_2 = np.matrix(((1,1,0,0),(1,1,1,0),(0,0,0,1),(0,0,0,0)))
    n_2 = np.array((1,1,2,1))
    phi_2 = configuration(S_2,n_2)

    #Mitochondrial data set
    n_mitochondrial = np.array((696,193,124,112,29,15,9,7,6,5,5,5,4,4,3,3,3,3,3,3,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1))
    S_mitochondrial = np.genfromtxt('S.csv',delimiter=',')
    S_mitochondrial = S_mitochondrial.astype(int)
    phi_mitochondrial = configuration(S_mitochondrial,n_mitochondrial)

    p_mitochondrial = np.zeros((3,5))
    p_1 = np.zeros((3,3))
    p_2 = np.zeros((3,3))

    b_1 = np.array((0,1,2))
    b_2 = np.array((3,4,5))
    b_mitochondrial = np.array((40,41,42,43,44))


    for i in range(3):
        for j in range(3):
            start_time = time.time()
            phiTable = configTable()
#            p_1[i][j] = prob(deepcopy(phi_1),theta[i],b_1[j],T,P)
#            p_2[i][j] = prob(deepcopy(phi_2),theta[i],b_2[j],T,P)
            p_mitochondrial[i][j] = prob(deepcopy(phi_mitochondrial),theta[i],b_mitochondrial[j],T,P)
            print (time.time()-start_time)
#    print p_1
#    print p_2
    print p_mitochondrial

#    start_time = time.time()
#    p_mitochondrial_base = prob(deepcopy(phi_mitochondrial),14,20,T,P)
#    print (time.time()-start_time)
#    print p_mitochondrial_base
