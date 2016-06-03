import thetaMLE
from csv2psi import csv2psi
import time
import sys
import numpy as np
from math import log10
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

if __name__ == '__main__':

    if len(sys.argv) < 3:

        print "Usage: python bemcmarking_thetaMLE.py  input_file_path output-file-path n_cores = 4 b = 1 ]"
        print " Will read a dataset from a .csv file, run thetaMLE (using the number of cores"
        print " given), and save a plot of the results to output-file-path."
        print "  e.g. \"python bemcmarking_thetaMLE.py  $FILE_PATH $OUTPUT_FILE_PATH 10 1\""
        print "  will run the code, with 1 mutation, on 10 cores."
    else:
        csv_path = str(sys.argv[1])
        out_path = str(sys.argv[2])
        try:
            S,nr,nc = csv2psi(csv_path)
        except Exception as e:
            print "Error occured reading from csv:",e

        if len(sys.argv) < 4:
            print 'No cores argument provided. Using the default number of cores for comparison: 4'
            cores_list = [4]
        else:
            cores = int(sys.argv[3])

        if len(sys.argv) < 5:
            print 'no b argument(s) provided. Uring standard value of 1'
            mutations_list = [1]
        else:
            mutations_list = map(int,sys.argv[4:])

        plt.figure(num = 1,figsize = (10,5), dpi = 100)

        color_list = ['black']
        style_list = ['-','--','-.']
        colors = {}
        styles = {}
        for i,mutations in enumerate(mutations_list):
            colors[mutations] = color_list[i%len(color_list)]
            styles[mutations] = style_list[i%len(style_list)]

        for mutations in mutations_list:
            print "Alalyzing dataset:\nS =\n%s\nnr =\n%s\nnc =\n%s\n Mutations = %i"%(str(S),str(nr),str(nc),mutations)

            if cores == 1:
                print 'Analyzing with 1 core...'
            else:
                print 'Analyzing with %i cores...'%cores

            t1 = time.time()
            mle,all_values = thetaMLE.thetaMLE(S, nr, nc, mutations = mutations, cores = cores)
            t2 = time.time()
            time_multiprocessing = t2 - t1

            all_values_flat = [(theta,likelihood) for s in all_values for (theta,likelihood) in s]
            all_values_flat.sort()

            print "Done (theta_MLE, log_theta_MLE = %.5f,%.5f). Elapsed time = %.2f sec.\n"%(mle,log10(mle),time_multiprocessing)
            #print "calculated (theta,L(theta))-values:"
            #print '\n'.join(['%f, %f'%(theta,likelihood) for theta,likelihood in all_values_flat])

            X_values = np.array([theta for theta,likelihood in all_values_flat])
            Y_values = np.array([likelihood for theta,likelihood in all_values_flat])

            plt.plot(X_values,Y_values,'|',color = colors[mutations], linestyle = styles[mutations],label = 'b = %i, $\\hat{\\theta}_{MLE} = %.8g$'%(mutations,mle))
            plt.axvline(mle,color = 'grey', linestyle='dotted')
        plt.xlabel(r'$\theta$')
        plt.ylabel(r'$q_{\theta , B \leq b}$')
        plt.legend()

        string_b_list = '-'.join(map(str,mutations_list))
        for suffix in ['.pdf','.eps','.jpg','.png']:
            full_path = '%s_b_%s_no_color%s'%(out_path,string_b_list,suffix)
            try:
                plt.savefig(full_path)
                print 'saved figure "%s"'%full_path
            except Exception as e:
                print "Error saving figure %s"%full_path,e
        plt.show()
        print "Shows over, folks!"
