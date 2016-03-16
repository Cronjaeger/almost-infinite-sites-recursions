import sys,time,re
import timeit
import cProfile, pstats, StringIO
import numpy as np
import math
from scipy.special import binom

import formula_outline as fo

#import simulation-code from coalescent-simulations
# CHANGE THE VALUE OF THE BELOW VARIABLE TO WHERE COALESCENT-SIMULATIONS
# HAS BEEN INSTALLED
coalescentSimulationsPath = '/home/mathias/programming/coalescent-simulations/'
#coalescentSimulationsPath = '/homes/cronjage/Programming/python/coalescent-simulations'
sys.path.insert(0,coalescentSimulationsPath)
import finiteSitesModell_investigations as fsmi
from psiFromCSV import psiFromCSV

def mylog10(x):
    if x > 0:
        return math.log(x,10)
    else:
        return -float('inf')



##############
# Import datasets to be analyzed below:
#

simDataPath = '/home/mathias/programming/coalescent-simulations/simData'
#simDataPath = '/homes/cronjage/Programming/python/coalescent-simulations/simData'
#Path for where simulated data is stored. Format should be as outlined in
# comments of psiFromCSV.py (in coalescent-simulations-module)


csvFileNames = [
# 'psi__N_2_theta_1pt0000_L_2.csv',\
# 'psi__N_2_theta_1pt0000_L_8.csv',\
# 'psi__N_2_theta_1pt0000_L_32.csv',\
# 'psi__N_8_theta_1pt0000_L_2.csv',\
# 'psi__N_8_theta_1pt0000_L_8.csv',\
# 'psi__N_8_theta_1pt0000_L_32.csv',\
# 'psi__N_16_theta_1pt0000_L_2.csv',\
# 'psi__N_16_theta_1pt0000_L_8.csv',\
# 'psi__N_16_theta_1pt0000_L_32.csv',\
# 'psi__N_32_theta_1pt0000_L_2.csv',\
# 'psi__N_32_theta_1pt0000_L_8.csv',\
# 'psi__N_32_theta_1pt0000_L_32.csv',\
# 'psi__N_2_theta_10pt0000_L_2.csv',\
# 'psi__N_2_theta_10pt0000_L_8.csv',\
# 'psi__N_2_theta_10pt0000_L_32.csv',\
# 'psi__N_8_theta_10pt0000_L_2.csv',\
# 'psi__N_8_theta_10pt0000_L_8.csv',\
# 'psi__N_8_theta_10pt0000_L_32.csv',\
# 'psi__N_16_theta_10pt0000_L_2.csv',\
# 'psi__N_16_theta_10pt0000_L_8.csv',\
# 'psi__N_16_theta_10pt0000_L_32.csv',\
# 'psi__N_32_theta_10pt0000_L_2.csv',\
# 'psi__N_32_theta_10pt0000_L_8.csv',\
# 'psi__N_32_theta_10pt0000_L_32.csv',\
# 'psi__N_32_theta_1pt0000_L_1024.csv',\
#
# 'hand-crafted/fromHeinsBook_L_4.csv',\
# 'hand-crafted/fromHeinsBook_L_10.csv',\
# 'hand-crafted/fromHeinsBook_L_1000.csv',\
# 'hand-crafted/fromHeinsBook_incompatibleColumn-pair_L_4.csv',\
# 'hand-crafted/fromHeinsBook_incompatibleColumn-pair_L_10.csv',\
# 'hand-crafted/fromHeinsBook_incompatibleColumn-pair_L_1000.csv',\
# 'hand-crafted/fromHeinsBook_with_a_2_L_4.csv',\
# 'hand-crafted/fromHeinsBook_with_a_2_L_10.csv',\
# 'hand-crafted/fromHeinsBook_with_a_2_L_1000.csv',\
# 'hand-crafted/hammer1995_0123_L_3.csv',\
# 'hand-crafted/hammer1995_0123_L_10.csv',\
# 'hand-crafted/hammer1995_0123_L_1000.csv',\
# 'hand-crafted/noMutations_N_5_L_4.csv',\
# 'hand-crafted/noMutations_N_5_L_10.csv',\
 'hand-crafted/noMutations_N_5_L_1000.csv',\
# 'hand-crafted/noMutations_N_5_L_1000000000.csv',\
 #'hand-crafted/aquardoAndGreenberg_fig1_0123_L_49.csv',\
# 'hand-crafted/aquardoAndGreenberg_fig1_0123_L_139500000.csv',\
 # 'hand-crafted/WardEtAl1991_0123_L_27.csv',\
]


#csvFilePaths = ['%s/%s'%(simDataPath,filename) for filename in csvFileNames]


p1 = re.compile('_L_\d+')
p2 = re.compile('\d+')
def getL(fileName):
    m1 = p1.search(fileName)
    m2 = p2.search(m1.group())
    return int(m2.group())

# def getTheta(filename)
#     return 1.0

#########
# Compute and output statistics of each dataset sequentially
#

for fileName in csvFileNames:

    filePath = '%s/%s'%(simDataPath,fileName)

    L = getL(fileName)

    theta = 10**-3 * L

    print "="*80
    print "Calculations for simulated dataset: %s"%fileName

    #read the dataset
    S,n = psiFromCSV(filePath,tidyUp = True)

    #print the dataset to standard output
    for i in xrange(S.shape[0]):
        print S[i],n[i]

    individuals = sum(n)
    haplotypes, segSites = S.shape
    non0Entries = sum(sum(S != 0))

    column0123count = [(sum(S[:,j] == 0), sum(S[:,j] == 1), sum(S[:,j] == 2),\
     sum(S[:,j] == 3) ) for j in range(S.shape[1])]
    non0typeCount = lambda t: sum(x != 0 for x in t[1:])

    oneTypeColumns = filter(
     lambda j: non0typeCount(column0123count[j]) == 1, range(S.shape[1]) )

    twoTypeColumns = filter(
     lambda j: non0typeCount(column0123count[j]) == 2, range(S.shape[1]) )

    threeTypeColumns = filter(\
     lambda j: non0typeCount(column0123count[j]) == 3, range(S.shape[1]) )

    #Compute list of incompatible columns-pairs#
    incompatibleColumns = fsmi.inconsistentColumnPairs([1]*S.shape[1],S)


    print "Basic statistics:"
    print " Individuals:             %i"%individuals
    print " Haplotypes:              %i"%haplotypes
    print " Segregating sites:       %i"%segSites
    # print "    Non-0 entires:           %i"%non0Entries
    if segSites > 1:
        print " Incompatible site-pairs: %i (%.0f%%)"%(len(incompatibleColumns), 200.0 * len(incompatibleColumns)/(segSites * (segSites - 1) ) )
    else:
        print" Incompatible site-pairs: 0"

    print '\n Excess mutation-statistics:'

    """
    Output count of number of types present in individual columns,
    e.g.
    ...1 non-0-type-columns count: 12
    ...2 non-0-type-columns count: 9
    ...3 non-0-type-columns count: 2
    """
    print '  1 non-0-type-columns count: %i\n  2 non-0-type-columns count: %i\n  3 non-0-type-columns count: %i\n'%(len(oneTypeColumns), len(twoTypeColumns), len(threeTypeColumns))

    """
    Output column-indices sorted by how many non-0-types they exhibit

    e.g.
    ...1 non-0-type columns: [0, 1, 2, 3, 4, 5, 7, 9, 13, 14, 15, 16]
    ...2 non-0-type columns: [8, 10, 11, 17, 18, 19, 20, 21, 22]
    ...3 non-0-type columns: [6, 12]
    """
    print '  1 non-0-type columns: %s\n  2 non-0-type columns: %s\n  3 non-0-type columns: %s\n'%(str(oneTypeColumns),str(twoTypeColumns),\
str(threeTypeColumns))

    """
    Output an incompatibility table, marking if two columns are incompatible
    e.g.
    ...    0  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20 21 22
    ... 0:                   x     x
    ... 1:                   x                                   x
    ... 2:
    ... 3:
    ... 4:                   x           x     x                 x  x        x
    ... 5:                   x     x        x  x              x     x        x
    ... 6:                         x     x  x  x           x  x  x  x  x  x  x
    ... 7:
    ... 8:                                        x           x  x           x
    ... 9:
    ...10:                                  x  x  x           x        x  x
    ...11:                                     x  x           x  x           x
    ...12:                                        x           x  x     x  x  x
    ...13:                                                       x        x
    ...14:
    ...15:
    ...16:                                                       x
    ...17:                                                       x  x     x
    ...18:                                                          x  x  x  x
    ...19:                                                             x  x  x
    ...20:                                                                   x
    ...21:                                                                   x
    ...22:
    """
    print 'Incompatibility-table:'
    incompTable = np.zeros((S.shape[1],S.shape[1]),dtype=bool)

    print '    '+' '.join(['%i '%i if i<10 else '%i'%i for i in range(S.shape[1])])
    for t in incompatibleColumns:
        incompTable[t[0],t[1]] = True
    for j1 in range(S.shape[1]):
        s = '%i:%s'%(j1,' '.join([" x" if x else "  " for x in incompTable[j1]]))
        if j1 < 10:
            s = ' ' + s

        print s
    #print an empty line after the table
    print ''



    print '- '+' - '*26
    print 'Benchmarks'
    print '  Presumptions:'
    print '      L     = %i'%L
    print '      theta = %f'%theta
    print '      P(i,j) = 1/3 * int(i!=j) '

    segSites = S.shape[1]
#    bPlus_list = [3,]
#    bPlus_list = [-1,0,1,2]
    bPlus_list = [0,1,2]
    b_list = [segSites + i for i in bPlus_list]


#    times = ()

    for b in b_list:

        print " b = %i: "%b

#        pr = cProfile.Profile()
#        pr.enable()

        t1 = time.time()
        prob,size = fo.prob_External(S,n,L,b,theta = theta)
#        prob,size = 0.05,5
        t2 = time.time()

        print "    probability:          %.20f"%prob
	print "    log10(probability) : %f"%mylog10(prob)
        print "    elapsed time:        %.3f sec."%(t2-t1)
        print "    config-table size:   %i"%size

#        s = StringIO.StringIO()
#        sortby = 'cumulative'
#        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
#        ps.print_stats()
#        # print " . "*10
#        print "\n Begining of profiler output: "
#        print s.getvalue()
#        print " End of profiler output"
#        print " . "*26
