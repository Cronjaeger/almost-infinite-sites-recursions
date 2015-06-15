# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 20:40:35 2015

@author: AleAviP, Mathias
"""

#from configurationTable import configTable,node
from configurations import configuration,coalesce,invisible,mutation
import numpy as np
import probabilities as proba
from copy import deepcopy

phi_null = configuration( S = np.matrix((0,)) , n = np.array((1,)) )

def prob(phi,theta,b,T,P):
    
    if ((np.all(phi.S==0)) and np.all((phi.n==1))):
        return 1
    
    coalescentIndicees = (i for i in range(len(phi.n)) if phi.n[i]>1)
    invisibleIndicees = ((i,k) for k in xrange(1,4) for i in xrange(phi.S.shape[0]))
    mutationIndicees = ((i,j,k) for i in xrange(phi.S.shape[0]) for j in xrange(phi.L) for k in xrange(4) if k != phi.S[i,j])

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
            if (phi.S[i,j]!=np.matrix((0,))):
                phi_mut = deepcopy(phi)
                fac_m = proba.p_mutation(phi_mut,i,j,k,theta,T,P)
                pro_m = prob( mutation(phi_mut,i,j,k), theta , b - 1, T, P )
                line_3 += fac_m * pro_m

    print ""
    print phi,b
    print line_1,fac_c,pro_c,len(list(coalescentIndicees))
    print line_2,fac_i,pro_i,len(list(invisibleIndicees))
    print line_3,fac_m,pro_m,len(list(mutationIndicees))
    p = line_1 + line_2 + line_3
#    print "womp"
    return p