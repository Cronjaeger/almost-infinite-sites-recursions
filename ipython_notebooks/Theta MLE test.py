
# coding: utf-8

# In[3]:

import numpy as np
import sys


# In[4]:

projectPath = '/home/mathias/programming/almost-infinite-sites-recursions'
sys.path.insert(0,projectPath)
from thetaMLE import thetaMLE


# In[5]:

S1 = np.matrix([0])
nR1 = np.array([2])
nC1 = np.array([1])
nC1_1 = np.array([1000])

S2 = np.matrix([[1,0],[0,0]])
nC2 = np.array([1,20])
nR2 = np.array([1,1])

S3 = np.matrix([
        [0, 0, 1, 0, 2, 0],
        [1, 2, 0, 0, 1, 0],
        [1, 1, 0, 1, 1, 0],
        [0, 0, 0 ,0, 1, 0]])
nC3 = np.array(
        [1, 1, 1, 1, 1, 9995])
nR3 = np.array([3,1,1,2])


# In[6]:

S3


# In[6]:

#MLE1 = thetaMLE(S1,nR1,nC1,1,verbose=False)
#MLE2 = thetaMLE(S2,nR2,nC2,2,verbose=True)
MLE3 = thetaMLE(S3,nR3,nC3,1,verbose=True)

