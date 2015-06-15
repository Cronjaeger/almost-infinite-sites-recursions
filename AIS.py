# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 16:27:45 2015

@author: AleAviP, Mathias
"""

from configurations import configuration,coalesce,invisible,mutation
import numpy as np
import probabilities as proba
import configurationTable as table

#TODO CHANGE tables 
#stateSpace = table.configTable()
#phi_knot=configuration(np.matrix(((0))),np.array([1]))
#stateSpace.add(phi_knot,p=1)

#WORKS!!!!!
#def AIS(phi,theta,b,T,P):
#    #Change initial conditions
#    probability=0
#    if ((np.all(phi.S==0)) and np.all((phi.n==1))):
#        return 1
#    else:
#        if (b>=0):
#            for i in (0,(np.shape(phi.S)[0]-1)):
#                if phi.n[i]>=2:
#                    phi_coal=phi
#                    probability = probability + proba.p_coalesce(phi_coal,i,theta)*AIS(coalesce(phi_coal,i),theta,b,T,P)
#                    return probability
#                else:
#                    pass
                
#TABLES DOESN'T work
#def AIS(phi,theta,b,T,P):
#    #Change initial conditions
#    probability=0
#    if ((np.all(phi.S==0)) and np.all((phi.n==1))):
#        return 1
#    else:
#        if (b>=0):
#            for i in (0,(np.shape(phi.S)[0]-1)):
#                if phi.n[i]>=2:
#                    phi_coal=phi
#                    if (stateSpace.contains(coalesce(phi_coal,i)) if coalesce(phi_coal,i).n[i]>1 else False):
#                        probability = probability+ proba.p_coalesce(phi_coal,i,theta)*stateSpace.get_p(coalesce(phi_coal,i))
#                    else:
#                        probability = probability + proba.p_coalesce(phi_coal,i,theta)*AIS(coalesce(phi_coal,i),theta,b,T,P)
#                        stateSpace.add(coalesce(phi_coal,i),p=probability)
#                    return probability
#                else:
#                    return probability

#def AIS(phi,theta,b,T,P):
#    #Change initial conditions
#    probability=0
#    if ((np.all(phi.S==0)) and np.all((phi.n==1))):
#        return 1
#    else:
#        if (b>=0):
#            for i in (0,(np.shape(phi.S)[0]-1)):
#                if phi.n[i]>=2:
#                    phi_new=coalesce(phi,i)
#                    probability = probability + proba.p_coalesce(phi,i,theta)*AIS(phi_new,theta,b,T,P)
#                    print probability
#                    return probability
#                else:
##                    pass
#                    if (b>0):
#                        stop= (b-phi.L-1) if ((b-phi.L)>0) else 0
#                        #Mutations + Invisible Mutations
#                        for j in (0,(phi.L+stop)):
#                            if (j<=(phi.L-1)):
#                                for k in (0,3):
#                                    if (phi.S[i,j]==k):
#                                        return 0.
#                                    else:
#                                        print phi
#                                        phi_mut = rmZeros(mutation(phi,i,j,k)) if (phi.L>1) else mutation(phi,i,j,k)
#                                        print phi_mut
#                                        probability = probability + proba.p_mutation(phi,i,j,k,theta,T,P)*AIS(phi_mut,theta,(b-1),T,P)
#                                        print probability
#                                        return probability
#                            else:
#                                return 0.
#                    else:
#                        return 0.

#Adding mutation
def AIS(phi,theta,b,T,P):
    #Change initial conditions
    probability=0
    if ((np.all(phi.S==0)) and np.all((phi.n==1))):
        return 1
    else:
        #COALESCE
        if (b>=0):
            for i in (0,(np.shape(phi.S)[0]-1)):
                if phi.n[i]>=2:
                    phi_coal=phi
                    probability = probability + proba.p_coalesce(phi_coal,i,theta)*AIS(coalesce(phi_coal,i),theta,b,T,P)
                else:
                    probability = probability
        #MUTATIONS and INVISIBLE states        
        elif (b>0):
            stop= (b-phi.L-1) if ((b-phi.L)>0) else 0
            #ALL configurations
            for i in (0,(np.shape(phi.S)[0]-1)):
                #ALL segregating sites
                for j in (0,(phi.L+stop)):
                #OLD segregating site                    
                    if (j<=phi.L):
                        #All different configurations
                        for k in (0,3):
                            if (phi.S[i,j]==k):
                                probability = probability
                            else:
                                phi_mut = phi
                                probability = probability + proba.p_mutation(phi_mut,i,j,k,theta,T,P)*AIS(mutation(phi_mut,i,j,k),theta,(b-1),T,P)
                    #NEW segregating site
                    else:
                        for k in (0,3):
                            phi_extra = phi
                            probability = probability + proba.p_invisible(phi_extra,i,k,theta,T,P)*AIS(invisible(phi_extra,i,k),theta,(b-1),T,P)
        else:
            probability = probability
    return probability
                            
#def runTest():
#    phi = configuration(np.matrix(((1,1),(1,0))),np.array((1,3),dtype = int))
#    p=float (1)/(3)
#    P=np.matrix(((0,p,p,p),(p,0,p,p),(p,p,0,p),(p,p,p,0)))
#    AIS(phi,2,2,100,P)
# 
#if __name__ == '__main__':
#    runTest()
    
                       
                        