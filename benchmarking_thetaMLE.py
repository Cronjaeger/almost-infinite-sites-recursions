import thetaMLE
import time
import sys
import numpy as np

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print 'No argument provided. Using the default number of cores for comparison: 4'
        cores = 4
    else:
        cores = int(sys.argv[1])
    print 'Benchmarking 1 vs. %i cores'%cores

    '''
    Our test-dataset is from Hammer 1995
    '''
    S = np.matrix([[0, 0, 3, 0],[0, 0, 0, 3],[0, 0, 0, 0],[3, 0, 0, 3]])
    nr = np.array([1, 1, 11, 3])
    nc = np.array([1, 997, 1, 1])
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
