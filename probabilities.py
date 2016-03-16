# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 16:14:02 2015

@author: AleAviP, Mathias
"""

#p=float (1)/(3)
#P=np.matrix(((0,p,p,p),(p,0,p,p),(p,p,0,p),(p,p,p,0)))

from configurations import configuration
import numpy as np

# class prob(configuration):
#
#     def hasRow(self,new_row):
#         for row in self.S:
#             if np.all(new_row == row):
#                 return True
#         return False
#
#     def whichRow(self,new_row):
#         for i,row in enumerate(self.S):
#             if np.all(new_row == row):
#                 return i
#         return 0

def p_coalesce(phi,i,theta):
    if phi.n[i]<2:
        return 0.
    else:
        return float (sum(phi.n)-1)/(sum(phi.n)-1+theta) * (phi.n[i]-1)/(sum(phi.n)-1)

def p_coalesce_new(phi,i,theta):
    if phi.nR[i]<2:
        return 0.
    else:
        return float(sum(phi.nR)-1)/(sum(phi.nR)-1+theta) * float(phi.nR[i]-1)/(sum(phi.nR)-1)


def p_mutation(phi,i,j,k,theta,T,P):

    #Compute new row
    new_S = np.matrix(phi.S)
    new_S[i,j]=k
    new_row = new_S[i,:]

    #Compare and compute S_new
         #EXISTENT
    if phi.hasRow(new_row):
        n_r=phi.n[phi.whichRow(new_row)]
        return float (theta)/(sum(phi.n)-1+theta) * P[k,phi.S[i,j]] * (phi.L)/(T) * (phi.n[i] * (n_r + 1))/(sum(phi.n))
    else:
        return float (theta)/(sum(phi.n)-1+theta) * P[k,phi.S[i,j]] * (phi.L)/(T) * (1)/(sum(phi.n))

def p_mutation_new(phi,theta,P,X,Y,Nr,Nc):
    N = sum(phi.nR)
    T = sum(phi.nC)
    p_mut_vs_coal = float(theta) / (N-1+theta)
    combinatorialFactor = float(Nr * Nc)/(N*T)
    transitionProb = P[X,Y]
    return p_mut_vs_coal * combinatorialFactor * transitionProb


def p_invisible(phi,i,k,theta,T,P):
    if (np.all(phi.S==0)):
        return float (theta)/(sum(phi.n)-1+theta) * P[k,0] * (phi.n[i])/(sum(phi.n))
    else:
        return float (theta)/(sum(phi.n)-1+theta) * P[k,0] * (T - phi.L)/(T) * (phi.n[i])/(sum(phi.n))


def runTests():
    print "It's working"

#runTests()
