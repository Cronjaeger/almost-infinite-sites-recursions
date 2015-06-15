# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 14:40:17 2015

@author: mathias & alejandra
"""

from configurationTable import configTable,node
from configurations import configuration,coalesce,invisible,mutation
import numpy as np
import probabilities as proba
from copy import deepcopy

phiTable = configTable()
phi_null = node( S = np.matrix((0,)) , n = np.array((1,)) )
#for b in xrange(10):
#    phiTable.add(phi_null,p=1.0,b=0)

def prob(phi,theta,b,T,P):

    if phiTable.contains(phi,b):
#        print "Wee!"
        return phiTable.get_p(phi,b)

    if phi == phi_null:
#        print "Woo"
        phiTable.add( phi = phi , p = 1.0 , b = b)
        return 1.0

    coalescentIndicees = (i for i in range(len(phi.n)) if phi.n[i]>1)
    invisibleIndicees = ((i,k) for k in xrange(1,4) for i in xrange(phi.S.shape[0]))
    mutationIndicees = ((i,j,k) for i in xrange(phi.S.shape[0]) for j in xrange(phi.L) for k in xrange(4) if k != phi.S[i,j] )

    fac_c,pro_c,fac_i,pro_i,fac_m,pro_m = 0,0,0,0,0,0


    line_1 = 0
    for i in coalescentIndicees:
        phi_coal = deepcopy(phi)
        fac_c = proba.p_coalesce(phi_coal,i,theta)
        pro_c = prob( coalesce(phi_coal,i), theta , b , T, P)
        line_1 += fac_c * pro_c

    line_2 = 0
    if b > 1:
        for i,k in invisibleIndicees:
            phi_inv = deepcopy(phi)
            fac_i = proba.p_invisible(phi_inv,i,k,theta,T,P)
            pro_i = prob( invisible(phi_inv,i,k) , theta , b - 1, T , P )
            line_2 += fac_i * pro_i

    line_3 = 0
    if b > 0:
        for i,j,k in mutationIndicees:
            if ((np.all(phi.S==0))==False):
                phi_mut = deepcopy(phi)
                fac_m = proba.p_mutation(phi_mut,i,j,k,theta,T,P)
                pro_m = prob( mutation(phi_mut,i,j,k), theta , b - 1, T, P )
                line_3 += fac_m * pro_m

#    print ""
#    print phi,b
#    print line_1,fac_c,pro_c,len(list(coalescentIndicees))
#    print line_2,fac_i,pro_i,len(list(invisibleIndicees))
#    print line_3,fac_m,pro_m,len(list(mutationIndicees))
    p = line_1 + line_2 + line_3
    phiTable.add( phi = phi , p = p , b = b)
#    print "womp"
    return p

#if __name__ == "__main__":
#
#    S_1 = np.matrix((0,))
#    n_1 = np.array((5,))
#    phi_1 = configuration(S_1,n_1)
#
#    P = (np.ones((4,4)) - np.eye(4))/3
#    theta = 4.0
#    T = 10
#
#    S_2 = np.matrix("0;1")
#    n_2 = np.array((2,1))
#    phi_2 = configuration(S_2,n_2)
##
##
##     p_1 = [prob(deepcopy(phi_1),theta,b,T,P) for b in range(5)]
##     phiTable = configTable()
##     p_2 = [prob(deepcopy(phi_2),theta,b,T,P) for b in range(5)]
##     print p_1
##     print phiTable.get_all_configurations()
##     print p_2
##
#    n_3 = np.array((1,1))
#    phi_3 = configuration(S_2,n_3)
#    phiTable = configTable()
#    p_3 = [prob(deepcopy(phi_3),theta,b,T,P) for b in range(5)]
#    print p_3
##
#    n_4 = np.array((2,))
#    phi_4 = configuration(S_1,n_4)
#    phiTable = configTable()
#    p_4 = [prob(deepcopy(phi_4),theta,b,T,P) for b in range(5)]
#    print p_4
#

#def runInconsistent():
#    S_i = np.matrix(((1,0),(1,0),(1,1)))
#    n_i = np.array((1,1,1))
#    phi_i = configuration(S_i,n_i)
#
#    theta = 2.0
#    T = 10
#    P = (np.ones((4,4)) - np.eye(4))/3
#
#    phiTable = configTable()
#    p_i = [prob(deepcopy(phi_i),theta,b,T,P) for b in range(5)]
#    print p_i
#
#runInconsistent()
#
#def typicalExample():
#    S_1 = np.matrix(((1,1,0,0),(1,1,1,0),(0,0,0,1),(0,0,0,0)))
#    n_1 = np.array((1,1,2,1))
#    phi_1 = configuration(S_1,n_1)
#
#    theta = 2.0
#    T = 10
#    P = (np.ones((4,4)) - np.eye(4))/3
#
#    phiTable = configTable()
#    p_1 = [prob(deepcopy(phi_1),theta,b,T,P) for b in xrange(3,6)]
#    print p_1
#
#typicalExample()
#
#def ChangingTheta():
#    list_theta = [i for i in range (101)]
#    S_2 = np.matrix("0;1")
#    n_2 = np.array((2,1))
#    phi_2 = configuration(S_2,n_2)
#
#    T=10
#    P = (np.ones((4,4)) - np.eye(4))/3
#    p_t1=[0 for j in range(101)]
#    for j in range(101):
#        phiTable = configTable()
#        p_t1[j] = prob(phi_2,float (list_theta[j])/(10),1,T,P)
#    print p_t1
#
#    p_t2=[0 for j in range(101)]
#    for j in range(101):
#        phiTable = configTable()
#        p_t2[j] = prob(phi_2,float (list_theta[j])/(10),2,T,P)
#
#    p_t3=[0 for j in range(101)]
#    for j in range(101):
#        phiTable = configTable()
#        p_t3[j] = prob(phi_2,float (list_theta[j])/(10),3,T,P)
#
#    T=100
#    p_t4=[0 for j in range(101)]
#    for j in range(101):
#        phiTable = configTable()
#        p_t4[j] = prob(phi_2,float (list_theta[j])/(10),1,T,P)
#
#    p_t5=[0 for j in range(101)]
#    for j in range(101):
#        phiTable = configTable()
#        p_t5[j] = prob(phi_2,float (list_theta[j])/(10),2,T,P)
#
#    p_t6=[0 for j in range(101)]
#    for j in range(101):
#        phiTable = configTable()
#        p_t6[j] = prob(phi_2,float (list_theta[j])/(10),3,T,P)
#
#    T=2
#    p_t7=[0 for j in range(101)]
#    for j in range(101):
#        phiTable = configTable()
#        p_t7[j] = prob(phi_2,float (list_theta[j])/(10),1,T,P)
#
#    p_t8=[0 for j in range(101)]
#    for j in range(101):
#        phiTable = configTable()
#        p_t8[j] = prob(phi_2,float (list_theta[j])/(10),2,T,P)
#
#    p_t9=[0 for j in range(101)]
#    for j in range(101):
#        phiTable = configTable()
#        p_t9[j] = prob(phi_2,float (list_theta[j])/(10),3,T,P)
#
#    p_i5 = prob(deepcopy(phi_i),theta,5,10,P)

    #Undone
S_i = np.matrix(((1,0),(1,0),(1,1)))
n_i = np.array((1,1,1))
phi_i = configuration(S_i,n_i)
theta = 2.0
T = 10
P = (np.ones((4,4)) - np.eye(4))/3

for b in range(1,8):
    phiTable = configTable()
    p_i = prob(deepcopy(phi_i),theta,b,100,P)
    print b,phiTable.get_size(),p_i

#    p_i3T100 = prob(deepcopy(phi_i),theta,3,100,P)
#
#    p_i3T1000 =prob(deepcopy(phi_i),theta,3,1000,P)
#
#    p_ib=[0 for j in range(101)]
#    for j in range(101):
#        phiTable = configTable()
#        p_ib[j] = prob(phi_i,float (list_theta[j])/(10),3,T,P)
#
#    baux=np.array((5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10))
#
#    for j in range(10):
#        phiTable = configTable()
#        p_ib[50+j*5] = prob(phi_i,baux[j],3,T,P)