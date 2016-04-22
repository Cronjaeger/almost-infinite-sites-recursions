import sys,os
import numpy as np
import time
import math
from benchmarking_singleExecution import waterson
from csv2psi import csv2psi
from formula_outline import prob_External

'''
This script will print a latex-table of use-cases intednded to help a
practitioner eyeball the computational time required to compute likelihoods for
simple tasks.
'''

def printTable(cesesDirectory,b_extras = [0]):
    '''
    Prints a table for each case (encoded as a .csv-file) in cesesDirectory.
    '''
    file_list = os.listdir(casesDirectory)
    csv_file_list = filter(lambda x: x[-4:] == '.csv',file_list)

    table_begin = r'\begin{table} \begin{tabular}{c|rrrrrrrrrrr}'
    table_end = r'\end{tabular} \end{table}'
    column_titles = [ 'Conf.',
                    #   'Sequences',
                    #   'Sites',
                    #   r'$b_{\min}$',
                      r'$b$',
                      r'$q_{\theta_w ; B \leq b}$',
                    #   r'$\log_10(q_{\theta_w}(\cdot , B \leq b))$',
                      r'$t_{b}$ (sec.)',
                      r'$c$',#'distinct configurations encountered',
                      r'$c_{>0}$',#'distinct configurations w. $q>0$',
                      r'$C$',#'terms total',
                      r'$C_{>0}$'#'number of non-0 terms'
                      ]

    first_row = vectorToRow(column_titles)

    table_string = '%s\n%s'%(table_begin,first_row)

    #b_extras = [b]

    for file_name in csv_file_list:
        file_path = '%s%s'%(casesDirectory,file_name)
        #print file_path
        S,nr,nc = csv2psi(file_path)

        n = sum(nr)
        L = sum(nc)

        deviants = set((1,2,3))
        b_min = sum( [ len(deviants.intersection(set([S[i,j] for i in xrange(S.shape[0])]))) * nc[j] for j in xrange(S.shape[1]) ] )

        theta_w = waterson(b_min,n)

        for b_extra in b_extras:
            b = b_min+b_extra
            #print S,nr,nc,'\n'
            t1 = time.time()
            prob,table = prob_External(S,nr,nc,b,theta_w)
            t2 = time.time()
            #runtime = '%.3f'%(t2 - t1)
            runtime = float(t2 - t1)
            size_small, size = table.get_size()
            size_smaller = sum([int( 0 < sum([1 for y in x.values() if y>0.0]) ) for x in table.get_all_configurations().values()])
            non_0_terms = sum([sum([1 for y in x.values() if y>0.0]) for x in table.get_all_configurations().values()])

            entries = [configurationToLaTex(S,nr,nc)]
            entries += map(
                           format_numbers,
                           [
                            # n,
                            # L,
                            # b_min,
                            b,
                            prob,
                            # myLog(prob),
                            runtime,
                            size_small,
                            size_smaller,
                            size,
                            non_0_terms
                           ]
                          )

            new_row = vectorToRow(entries)

            print '\\\\\n%s'%new_row

            table_string += '\\\\\n%s'%new_row

    table_string += '\n%s'%table_end

    print table_string

def format_numbers(x):
    if type(x) == int:
        return str(x)
    else:
        return '$%.3g$'%x

def myLog(x):
    if x > 0:
        return math.log10(x)
    else:
        return float('-inf')

def configurationToLaTex(S,nr,nc):
    strings = map(matrixToBmatrix,[S,nr,nc])
    return '$\\left( %s , %s , %s \\right)$'%(matrixToBmatrix(S),matrixToBmatrix(nr),matrixToBmatrix(nc))

def vectorToBMatrix(vec):
    return '\\begin{bmatrix} %s \\end{bmatrix}'%r' \\ '.join([str(v) for v in vec])

def vectorToRow(vec):
    return r' & '.join([str(v) for v in vec])

def matrixToBmatrix(M):
    if M.ndim == 1:
        return vectorToBMatrix(M)
    else:
        row_strings = map(vectorToRow,M)
        return vectorToBMatrix(row_strings)

if __name__ == '__main__':
    if len(sys.argv) == 1:
        casesDirectory = str(os.getcwd())+'/'
    else:
        casesDirectory = str(sys.argv[1])

    if len(sys.argv) > 2:
        b_extras = map(int,sys.argv[2:])
    else:
        b_extras = [0]

    printTable(casesDirectory,b_extras)
