import formula_outline as fo
from csv2psi import csv2psi
from thetaMLE import waterson
import time
import sys
import numpy as np

def main():
    filename = sys.argv[1]
    extra_mutations = int(sys.argv[2])

    S,nR,nC = csv2psi(filename)
    print "Alalyzing dataset from \'%s\':\n\nS =\n%s\nnr =\n%s\nnc =\n%s\n"%(filename,str(S),str(nR),str(nC))
    #print S,Nr,Nc

    deviants = set((1,2,3))
    min_mutations = sum( [ len(deviants.intersection(set([S[i,j] for i in xrange(S.shape[0])]))) * nC[j] for j in xrange(S.shape[1]) ] )
    n = sum(nR)
    B_max = min_mutations + extra_mutations

    thetaGuess = waterson(min_mutations,n)

    print 'Weighted count of segregating sites = %i'%min_mutations
    print 'Maximal number of mutations in analysis = %i'%B_max
    print 'theta (waterson estimate) = %.3f'%thetaGuess

    print "running analysis ...",
    t1 = time.time()
    prob,table = fo.prob_External(S,nR,nC,B_max,thetaGuess)
    t2 = time.time()

    size = table.get_size()
    size_smaller = sum([int( 0 < sum([1 for y in x.values() if y>0.0]) ) for x in table.get_all_configurations().values()])
    non_0_terms = sum([sum([1 for y in x.values() if y>0.0]) for x in table.get_all_configurations().values()])
    t_diff = t2 - t1

    print 'done!'
    print 'q(S,nR,nC,B<%i) = %f '%(B_max+1,prob)
    print 'time = %.3f s'%t_diff
    print 'Encountered configurations = %i'%size[0]
    print 'Encountered configurations w pos. probability =%i'%size_smaller
    print 'Total number of terms considered = %i'%size[1]
    print 'Total number of positive terms = %i'%non_0_terms
    #print '(table_size,teble_size_extended) = %s'%str(size)
    #print str(table.get_all_configurations())

if __name__ == "__main__":
    main()
