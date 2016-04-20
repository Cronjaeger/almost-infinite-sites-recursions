import thetaMLE
from csv2psi import csv2psi
import time
import sys
import numpy as np
from math import log10

if __name__ == '__main__':

    if len(sys.argv) < 2:

        print "Usage: python bemcmarking_thetaMLE.py  file_path  [*n_cores = 4]"
        print " Will read a dataset from a .csv file and run thetaMLE once for each number of cores given."
        print "  e.g. \"python bemcmarking_thetaMLE.py  $FILE_PATH  1 10\" will run the code twice,"
        print "  - once with 1 core,"
        print "  - once with 10 cores."

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
            cores_list = [4]
        else:
            cores_list = [int(i) for i in sys.argv[2:]]
        #print 'Benchmarking 1 vs. %i cores'%cores


        # S  = np.matrix([[0,0],[0,1]])
        # nr = np.array([1 , 1])
        # nc = np.array([1 , 1])

        extra_mutations = 1

        print "Alalyzing dataset:\nS =\n%s\nnr =\n%s\nnc =\n%s\nExtra Mutations = %i"%(str(S),str(nr),str(nc),extra_mutations)

        # print "Alalyzing with 1 core..."
        # t1 = time.time()
        # mle1 = thetaMLE.thetaMLE(S, nr, nc, extra_mutations_allowed = extra_mutations, cores = 1)
        # t2 = time.time()
        # time_1 = t2 - t1
        # print mle1
        # print "Done (theta = %f). Elapsed time = %.2f sec.\n"%(mle1,time_1)

        for cores in cores_list:
            if cores == 1:
                print 'Analyzing with 1 core...'
            else:
                print 'Analyzing with %i cores...'%cores
            t1 = time.time()
            mle,all_values = thetaMLE.thetaMLE(S, nr, nc, extra_mutations_allowed = extra_mutations, cores = cores)
            t2 = time.time()
            time_multiprocessing = t2 - t1

            all_values_flat = [(theta,likelihood) for s in all_values for (theta,likelihood) in s]
            all_values_flat.sort()

            print "Done (theta,log_theta = %.5f,%.5f). Elapsed time = %.2f sec.\n"%(mle,log10(mle),time_multiprocessing)
            print "calculated (theta,L(theta))-values:"
            #print all_values
            print '\n'.join(['%f, %f'%(theta,likelihood) for theta,likelihood in all_values_flat])

        #print 'relative speedup (time 1 core / time %i cores) %f'%(cores,time_1/time_multiprocessing)
#
# def countMinMutationsWithoutRegardForIncompatiabilities(S,nc):
#     deviants = set((1,2,3))
#     minMutations = sum( [ len(deviants.intersection(set([S[i,j] for i in xrange(S.shape[0])]))) * nc[j] for j in xrange(S.shape[1]) ] )
#     return minMutations
