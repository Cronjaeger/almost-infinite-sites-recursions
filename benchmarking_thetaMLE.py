import thetaMLE
from csv2psi import csv2psi
import time
import sys
import numpy as np

if __name__ == '__main__':

    if len(sys.argv) < 2:

        print "Usage: python bemcmarking_thetaMLE.py  file_path  [n_cores = 4]"
        print " Will read a dataset from a .csv file and run thetaMLE twice."
        print "  - once with 1 core"
        print "  - once with n_cores cores"

        # '''
        # Our test-dataset is from Hammer 1995
        # '''
        # S = np.matrix([[0, 0, 3, 0],[0, 0, 0, 3],[0, 0, 0, 0],[3, 0, 0, 3]])
        # nr = np.array([1, 1, 11, 3])
        # nc = np.array([1, 997, 1, 1])
        #
        # print 'No path to test-dataset provided. test-data from Hammer et. al. (1995)'

        # print 'No path to test-dataset provided. Using the following:'
        # print 'S = %s'%str(s)
        # print 'nr = %s'%str(nr)
        # print 'nc = %s'%str(nc)

    else:
        csv_path = str(sys.argv[1])
        try:
            S,nr,nc = csv2psi(csv_path)
        except Exception as e:
            print "Error occured reading from csv:",e

        if len(sys.argv) < 3:
            print 'No seccond argument provided. Using the default number of cores for comparison: 4'
            cores = 4
        else:
            cores = int(sys.argv[2])
        print 'Benchmarking 1 vs. %i cores'%cores


        # S  = np.matrix([[0,0],[0,1]])
        # nr = np.array([1 , 1])
        # nc = np.array([1 , 1])

        extra_mutations = 2

        print "Alalyzing dataset:\nS =\n%s\nnr =\n%s\nnc =\n%s\nExtra Mutations = %i"%(str(S),str(nr),str(nc),extra_mutations)

        print "Alalyzing with 1 core..."
        t1 = time.time()
        mle1 = thetaMLE.thetaMLE(S, nr, nc, extra_mutations_allowed = extra_mutations, cores = 1)
        t2 = time.time()
        time_1 = t2 - t1
        print mle1
        print "Done (theta = %f). Elapsed time = %.2f sec.\n"%(mle1,time_1)

        print 'Analyzing with %i cores...'%cores
        t1 = time.time()
        mle2 = thetaMLE.thetaMLE(S, nr, nc, extra_mutations_allowed = extra_mutations, cores = cores)
        t2 = time.time()
        time_multiprocessing = t2 - t1
        print "Done (theta = %f). Elapsed time = %.2f sec.\n"%(mle2,time_multiprocessing)

        print 'relative speedup (time 1 core / time %i cores) %f'%(cores,time_1/time_multiprocessing)
#
# def countMinMutationsWithoutRegardForIncompatiabilities(S,nc):
#     deviants = set((1,2,3))
#     minMutations = sum( [ len(deviants.intersection(set([S[i,j] for i in xrange(S.shape[0])]))) * nc[j] for j in xrange(S.shape[1]) ] )
#     return minMutations
