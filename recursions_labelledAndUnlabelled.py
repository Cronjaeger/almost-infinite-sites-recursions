iimport numpy as np
from configurations import mutation_new, configuration_new
#from configurationTable import node as UUconfData

class Configuration(object):
    """
    A parent-class defining common fields and methods for configurations
    """
    def __init__(self, conf = None, fullSeqArr = None):
        '''
        Construct a configuration Must be passed either a full 2d array (rows =
        individuals, columns = aligned sites) or another Configuration object of
        the same type (in which case the constructor just copies)
        '''
        if conf is None:

            tempFullSeqArr = np.array(fullSeqArr, ndmin = 2)

            #compute the number of sequences
            n = tempFullSeqArr.shape[0]
            assert n > 0
            self.n = n

            #compute the number of sites
            L = tempFullSeqArr.shape[1]
            assert L > 0
            self.L = L

            # convert the array to whatever format encodes sequences in this
            # particular instance of the class.
            # Note: the method __seqs2config__(self) will raise a
            #       NotImplementedError by default.
            self.data = self.__seqs2config__(tempFullSeqArr)

        else:
            #we have been passes a configuration; the constructor should just
            # copy
            assert type(self) == type(conf)
            self.data = conf.data


    def __seqs2config__(self, fullSeqArr):
        raise NotImplementedError()

    def getSegSites(self):
        return NotImplemented

    def __hash__(self):
        # NOTE: If self.data is of an unhashable type, this method must be
        # re-implemented!
        return NotImplemented

    def __eq__(self,other):
        return NotImplemented


class LL_configuration(Configuration):
    '''
    A configuration with labelled sequences and labelled sites
    '''

    def __seqs2config__(self, fullSeqArr):
        return np.array(fullSeqArr , ndmin = 2)

    def __hash__(self):
        # converts the matrix self.data to a tuple of row-tuples, and hash the
        # resulting tuple of tuples
        return hash(tuple(map(tuple,self.data)))

    def __eq__(self,other):
        try:
            return np.all(self.data == other.data)
        except AttributeError:
            raise NotImplementedError("type(self) = %s, type(other) = %s"%( str(type(self)), str(type(other)) ) )

    # def __segsites__(self):
    #     typecountsByColumn =
    #       np.apply_along_axis( arr    = self.data,
    #                            axis   = 0,
    #                            func1d = lambda column: len(set(column)) )
    #     return ( x > 1 for x in typecountsByColumn )


class UL_configuration(Configuration):
    '''
    A configuration with unlabelled sequences and labelled sites
    '''
    def __seqs2config__(self, fullSeqArr):
        rowCounts = dict()
        rows = map(tuple, fullSeqArr)
        for row in rows:
            try:
                rowCounts[row] += 1
            except KeyError:
                rowCounts[row] = 1

        assert sum(rowCounts.values()) == self.n
        return rowCounts

    def __hash__(self):
        return hash( sorted(self.data.items()) )

    def __eq__(self,other):
        try:
            return sorted(self.data.items()) == sorted(other.data.items())
        except AttributeError:
            raise NotImplementedError("type(self) = %s, type(other) = %s"%( str(type(self)), str(type(other)) ) )

class LU_configuration(Configuration):
    '''
    A configuration with labelled sequences and unlabelled sites
    '''
    def __seqs2config__(self, fullSeqArr):
        colCounts = dict()
        columns = map(tuple, np.transpose(fullSeqArr))
        for col in columns:
            try:
                colCounts[col] += 1
            except KeyError:
                colCounts[col] = 1

        assert sum(colCounts.values()) == self.L
        return colCounts

    def __hash__(self):
        return hash( sorted(self.data.items()) )

    def __eq__(self,other):
        try:
            return sorted(self.data.items()) == sorted(other.data.items())
        except AttributeError:
            raise NotImplementedError("type(self) = %s, type(other) = %s"%( str(type(self)), str(type(other)) ) )

class UU_Configurations(Configuration):
    '''
    A configuration with unlabelled sequences and unlabelled sites.
    This is a wrapper around the class configuration_new from configurations.py
    '''
    def __seqs2config__(self, fullSeqArr):
        #compute S,nr,nc, and return configuration_new(S,nr,nc)

        rowCounts = dict()
        rows = map(tuple, fullSeqArr)
        for row in rows:
            try:
                rowCounts[row] += 1
            except KeyError:
                rowCounts[row] = 1

        uniqueRows = rowCounts.keys()
        nr = np.array(rowCounts.values())
        seqArr_noDuplicateRows = np.array(uniqueRows, ndmin = 2)

        colCounts = dict()
        columns = map(tuple, np.transpose(seqArr_noDuplicateRows))
        for col in columns:
            try:
                colCounts[col] += 1
            except KeyError:
                colCounts[col] = 1

        assert sum(colCounts.values()) == self.L

        uniqueCols = colCounts.keys()
        S = np.transpose( np.array(uniqueCols, ndmin = 2) )
        nc = np.array(colCounts.values())

        assert sum(nr) == self.n
        assert sum(nc) == self.L
        assert len(nr) == S.shape[0]
        assert len(nc) == S.shape[1]

        return configuration_new(S,nr,nc)

    def __hash__(self):
        return hash(self.data)

    def __eq__(self, other):
        try:
            return self.data == other.data
        except AttributeError:
            raise NotImplementedError("type(self) = %s, type(other) = %s"%( str(type(self)), str(type(other)) ) )

class AncestralGraphNode(object):

    def __init__(self, conf, parents, probabilities):
        assert hasattr(conf,'data')
        assert len(parents) == len(probabilities)
        assert all( p <= 1.0 for p in probabilities)
        self.conf = conf
        self.parents = parents
        self.probabilities = probabilities

    def __hash__(self):
        return hash(self.conf)

    def __eq__(self,other):
        if hasattr(conf,'data'):
            return self.conf == other.conf
        else:
            raise NotImplementedError("type(self) = %s, type(other) = %s"%( str(type(self)), str(type(other)) ) )

class AncestralGraph(object):
    '''
    A directed graph, where each node is a configuration and each arc connects
    configurations separated by a single event. Node-weights are likelihoods and
    arc-weights are (forward) trnasition probabilities.
    '''
    def __init__(self, model):
        '''
        Initialize an empty ancestral graph
        '''
        self.nodes2probs = {}
        self.arcss2probs  = {}

    def addNode(self, conf, prob = None):

        try:
            assert conf not in self.nodes2probs
        except AssertionError:
            raise LookupError('Configuration already exists')

        assert ( prob is None ) or ( prob >= 0.0 and prob <= 1.0 )
        conf.nodes2probs[conf] = prob

    def addArc(self, origin, terminus, prob):
        if not (prob >= 0.0 and prob <= 1.0):
            raise ValueError('arcs must be assigned a probability between 0 and 1')
        self.arcss2probs[(origin,terminus)] = prob

    def resolveLikelihood(self, conf):
        pass

    def getStaes(self):
        pass

class Model(object):
    '''
    An instance of a finite sites model with bounded number of mutations.
    '''
    def __init__(self,
                 theta,
                 substitution_probs,
                 boundary_condition,
                 mutation_budget):

        # self.__validateInput__()
        self.theta = theta
        self.pi = boundary_condition
        self.P = substitution_probs
        self.mutation_bound = mutation_budget
        self.n_nucleotides = P.shape[0]

    def __validateInput__(self):
        # Fill out if nessecary
        pass


    def iterM(self,conf):
        '''
        A generator iterating over M(conf)
        returns an iterator over all tuples (conf2, prob) such that:
         - conf is obtainable by mutating a single nucleotide in conf2
         - prob is the (forward) probability of this mutation ocurring next
            given conf2
           * note that we do rule out conf2 == conf1 here, to do that, make sure
             that the trnasition matrix of the model has no diagonals.
        '''
        n,L,K = conf.n, conf.L, self.n_nucleotides
        p_mut_base = self.theta / ( (n - 1 + self.theta) * n * L)
        if isinstance(conf, LL_configuration):
            #p_mut_base =  p_mut_base / (n * L)
            tempSeqData = np.array(conf.data)
            i, j, k = 0,0,0
            while i < n:
                while j < L:
                    x = conf.data[i,j]
                    while k < K:
                        if k != conf.data[i,j] and self.P[k,x] > 0:
                            tempSeqData[i,j] = k
                            newConf = LL_configuration(conf) # copy old configuration
                            newConf.data = np.array(tempSeqData) #update data
                            trans_prob = p_mut_base * self.P[k,x]
                            yield (newConf, trans_prob)
                        k += 1
                    tempSeqData[i,j] = x
                    j += 1
                i += 1

        elif isinstance(conf, UL_configuration):
            #tempRowData  = dict(conf.data)
            rowCounts = conf.data
            for row,n_row in rowCounts.iteritems():

                #generate all row-counts that can already be determined
                if n_row == 1: # do not include old row if it is a singleton-row
                    rowCountsNew_generator = tuple( (r,n) for r,n in rowCounts.iteritems() if r != row )
                else: # reduce row-count by one if the row is a not a singleton row
                    rowCountsNew_generator = tuple( (r,n - int(r==row)) for r,n in rowCounts.iteritems() )

                j,k = 0,0:
                while j < L:
                    x = row[j]
                    while k < K:
                        if self.P[k,x] > 0:
                            newRow = tuple( (k if j2 == j else s) for j2,s in enumerate(row) )
                            rowCountsNew = dict(rowCountsNew_generator)
                            try:
                                rowCountsNew[newRow] += 1
                            except KeyError:
                                rowCountsNew[newRow] = 1
                            trans_prob = p_mut_base * rowCountsNew[newRow] * self.P[k,x]
                            newConf = UL_configuration(conf = conf) #copy old configutation
                            newConf.data = rowCountsNew # update data
                            yield (newConf, p_trans)
                        k += 1
                    j += 1

        elif isinstance(conf, LU_configuration):
            #tempColData  = dict(conf.data)
            colCounts = conf.data
            for col,n_col in colCounts.iteritems():

                #generate all column-counts that can already be determined
                if n_col == 1: # do not include old column if it is a singleton-column
                    colCountsNew_generator = tuple( (c,n) for c,n in colCounts.iteritems() if c != col )
                else: # reduce column-count by one if the column is a not a singleton column
                    colCountsNew_generator = tuple( (c,n - int(c==col)) for r,n in colCounts.iteritems() )

                i,k = 0,0:
                while i < n:
                    x = col[i]
                    while k < K:
                        if self.P[k,x] > 0:
                            newCol = tuple( (k if i2 == i else s) for i2,s in enumerate(col) )
                            colCountsNew = dict(colCountsNew_generator)
                            try:
                                colCountsNew[newCol] += 1
                            except KeyError:
                                colCountsNew[newCol] = 1
                            newConf = LU_configuration(conf = conf)
                            newConf.data = colCountsNew
                            trans_prob = p_mut_base * colCountsNew[newCol] * self.P[k,x]
                            yield (newConf, trans_prob)
                        k += 1
                    i += 1

        elif isinstance(conf, UU_configuration):
            if model.n_nucleotides > 4:
                raise NotImplementedError('Implementation currently does not support models with >4 nucleotides.')
            phi = conf.data
            i,j,k = 0,0,0
            I,J = phi.S.shape[0],phi.S.shape[1]
            while i < I:
                while j < J:
                    while k < K:
                        if self.P[k,x] > 0:
                            phi_mut,Nr,Nc = mutation_new(phi,i,j,k)
                            newConf = LU_configuration(conf = conf)
                            newConf.data = phi_mut
                            trans_prob = p_mut_base * Nr, Nc * self.P[k,x]
                            yield (newConf, trans_prob)
                        k += 1
                    j += 1
                i+= 1

        else: #Having exhausted all known cases for the type of conf, raise an error.
            raise NotImplementedError('Iterating over mutation-ancestors is not implemented for configurations of type = %s'%str(type(conf)) )

    def iterC(self, conf):
        '''
        A generator iterating over C(conf)
        returns an iterator over all tuples (conf2, prob) such that:
         - conf2 is obtainable by merging two lineages of conf
         - prob is the (forward) probability of conf2 turning into conf as a
           result of the next event.
        '''
        p_coal_base = 1 / (theta + conf.n - 1)
        if isinstance(conf, LL_configuration):
            p_coal_base /= conf.n
            i1 = 0
            while i1 < n - 1:
                row1 = conf.data[i1,:]
                matchingRowIndices = tuple(i2 for i2 in xrange(i1 + 1, conf.n + 1) if np.all(row1 == conf.data[i2,:]) )
                if len(matchingRowIndices) > 0:
                    newConf1 = LL_configuration(conf = conf)
                    newConf1.n -= 1
                    newConf1.data = np.delete(arr = newConf1.data, obj = i1, axis = 0)
                for i2 in matchingRowIndices:
                    newConf2 = LL_configuration(conf = conf)
                    newConf2.n -= 1
                    newConf2.data = np.delete(arr = newConf2.data, obj = i2, axis = 0)
                    yield (newConf1, p_coal_base)
                    yield (newConf2, p_coal_base)
                i1 += 1

        elif isinstance(conf, UL_configuration):
            mergeableRows = filter(lambda x: x[1] > 1, conf.data.iteritems() )
            for row,nRow in mergeableRows:
                newConf = UL_configuration(conf = conf)
                newRowDict = dict(conf.data)
                newConf.n -= n - 1

                newRowdict.data[row] -= 1
                newConf.data = newRowDict

                transProb = p_coal_base * (nRow - 1) / float(n - 1)
                yield (newConf, transProb)

        elif isinstance(conf, LU_configuration):
            # TODO: READ THROUGH THIS AGAIN!
            p_transition = p_coal_base / conf.n
            seqMatrix_transposed = np.array( conf.data.keys() )
            assert seqMatrix.shape[1] == n
            colCounts = conf.data.values()
            while i1 < n - 1:
                matchingRowIndices = tuple(i2 for i2 in xrange(i1+1, n+1) if np.all(seqMatrix_transposed[ :, i1] == seqMatrix_transposed[:, i2]) )
                if len(matchingRowIndices) > 0:
                    conf1 = LU_configuration(conf = conf)
                    conf1.n -= 1
                    columns1 = (tuple( np.delete(column, i1) ) for column in seqMatrix_transposed )
                    conf1.data = dict(zip( columns1 , colCounts ) )
                    for i2 in matchingRowIndices
                        conf2 = LU_configuration(conf = conf)
                        conf2.n -= 1
                        columns2 = (tuple( np.delete(column, i2) ) for column in seqMatrix_transposed )
                        conf2.data = dict(zip( columns2, colCounts ))
                        yield (conf1, p_coal_base )
                        yield (conf2, p_coal_base )
                i1 += 1

        elif isinstance(conf, UU_configuration):
            pass # TODO: fill in
        else: #Having exhausted all known cases for the type of conf, raise an error.
            raise NotImplementedError('Iterating over coalecence-ancestors is not implemented for configurations of type = %s'%str(type(conf)) )

def recursive_likelihood(conf, model):

    n = conf.n
    theta = model.theta
    pi = model.pi
    b_extra = model.extra_mutations7


    #Base Case: conf contains a single sequence
    if n == 1:
        return pi(conf)


    #Recursion case: conf doesw not contain a single sequence

    C = model.generateCoalescenceParents(conf)
    coalescence_term = 0.0
    for x in C:
        x_prob = recursive_likelihood(x['conf'], model)
        transition_prob = x['prob']
        coalescence_term +=  x_prob * transition_prob

    M = model.generateMutationParents(conf)
    mutation_term = 0.0
    for x in M:
        x_prob = recursive_likelihood(x['conf'], model)
        transition_prob = x['prob']
        mutation_term +=  x_prob * transition_prob

    return coalescence_term + mutation_term
