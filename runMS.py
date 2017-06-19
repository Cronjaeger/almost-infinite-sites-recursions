''' A wrapper function for running MS, as well as auxiliary scripts for parsing
the output. '''

import subprocess
import numpy as np
import time

def run_sim(args):
    p = subprocess.Popen(args,shell = False, stdout=subprocess.PIPE)
    #p = subprocess.Popen(args,shell = False, stdout=subprocess.PIPE, stdin = subprocess.PIPE, )
    output = p.communicate()[0]
    #print 'run_sim_output: ',output
    return output

def generate_seed():
    '''
    an auxiliary function which allows us to generate a random seed to opass to
    ms, but which depends on the state of the ransom number generator used by
    np.random. This way, seeding the random number generator at the beginning of
    a script with "np.random.set_state()" will also guarantee that the seeds
    used in any sub-processes executring ms will be the same.
    '''
    return int(np.random.sample(1)*(2**31))

def ms_sim_theta(n, theta, N = 1, rho = 0.0 , nsites = 1000000, seed = generate_seed(), with_trees = True, extra_args = []):
    if n<=1 or theta <0 or N<1: raise(ValueError('invalid input to ms_sim_theta: n=%i, theta=%f, N=%i.'%(n,theta,N)))
    if with_trees:
        extra_args.append('-T')
    args = ['ms',str(n),str(N),'-seeds',str(seed),'-t',str(theta),'-rho', str(rho), str(nsites)]+extra_args
    output_raw = run_sim(args)
    output_parsed = parse_ms_output(output_raw)
    return output_parsed

def ms_sim_sites(n, s, N = 1, rho = 0.0, nsites = 1000000, seed = generate_seed(), with_trees = True, extra_args = []):
    if n<=1 or s <1 or N<1: raise(ValueError('invalid input to ms_sim_sites: n=%i, s=%i, N=%i.'%(n,s,N)))
    if with_trees:
        extra_args.append('-T')
    args = ['ms',str(n),str(N),'-seeds',str(seed),'-s',str(s),'-rho', str(rho), str(nsites)] + extra_args
    output_raw = run_sim(args)
    output_parsed = parse_ms_output(output_raw)
    return output_parsed

def parse_ms_output(raw_input):
    'auxiliary function for parsing raw output of ms'
    raw = str(raw_input)
    lines = str.split(raw,'\n')
    chunks =  str.split(raw,r'//')
    # print 'parsing raw input:'
    # print raw
    # print '='*80
    '''
    for refenrence, this is what a chunk looks like:
        ['',
         'segsites: 2',
         'positions: 0.3134 0.5345 ',
         '00',
         '00',
         '01',
         '11',
         '00',
         '',
         '']
    '''
    S_list = []
    trees_list = []
    positions_list = []

    for i,chunk in enumerate(chunks[1:]): #the first entry lists input,seed and trees (if applicable)

        # add an extra empty line to the last line (this way all chunks terminate in two empty lines)
        if chunk == chunks[-1]:
            chunk += '\n'

        #separate the chunk into lines
        chunk_lines = str.split(chunk,'\n')

        #the 0th line is empty

        #the next lines represent trees (if they are simulated)
        if chunk_lines[1][0] == 's':
            n_trees = 0
            trees = []
            #print 'NO trees. chunk_lines[1]=%s'%chunk_lines[1]
        else:
            n_trees = 1
            trees = [parse_tree_from_ms(chunk_lines[n_trees])]
            while chunk_lines[n_trees+1][0] != 's':
                n_trees += 1
                trees.append(parse_tree_from_ms(chunk_lines[n_trees]))

        trees_list.append(trees)

        #the following line contains the number of segregating sites
        segsites_str = chunk_lines[1+n_trees]
        n_segsites = parse_segsites(segsites_str)

        # when there are 0 segregating sites, there is nothing more to parse.
        if n_segsites == 0:
            positions = []
            S = np.array([0],dtype=int,ndmin=2)

        #when there are segregating sites, there is also output to parse
        else:
            positions_str = chunk_lines[2+n_trees]
            positions = parse_positions(positions_str)

            rows_str = chunk_lines[(3+n_trees):-2] # these are the lines which correspond to output
            rows_arr = map(parse_row,rows_str) #each row is converted to a 1d array of integers
            S = np.r_[rows_arr]

        positions_list += [positions]
        S_list.append(S)

    return {'raw':raw,
            'lines':lines,
            'metadata':parse_metadata(chunks[0]),
            'S_list':S_list,
            'positions':positions_list,
            'trees_list':trees_list,
            'n_segsites':n_segsites,
            'experiments':zip(S_list,positions_list,trees_list)}

def parse_segsites(input_str):
    # print 'parsing segsites from string: %s'%input_str
    return int(str.split(input_str)[1])

def parse_positions(input_str):
    # print 'parsing positions from string: %s'%input_str
    return map(float,str.split(input_str)[1:])

def parse_tree_from_ms(input_str):
    '''
    takes a tree-string like '[10](1:0.201,(2:0.078,3:0.078):0.123);', and
    outputs
    '''
    # print 'parsing tree from string: %s'%input_str
    if input_str[0]=='[':
        closing_bracket_position = input_str.index(']')
        sites = int(input_str[1:closing_bracket_position])
        newick_str = input_str[closing_bracket_position+1:]
    elif input_str[0]=='(':
        newick_str = input_str
        sites = 'all'
    else:
        raise TypeError("cannot parse a tree from input string '%s'"%input_str)

    return {'sites':sites,'str':newick_str}

def parse_row(row_str):
    return np.array([int(i) for i in row_str], dtype = int)

def parse_metadata(input_str):
    lines = str.split(input_str, '\n')
    return {'command':lines[0],'seed':int(lines[1])}
