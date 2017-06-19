'''implements the approach to counting histories outlined in my 3rd note,
available at https://www.overleaf.com/read/gbqwxpqhjkjy'''

from collections import Counter
from itertools import product, chain, combinations
import numpy as np
from runMS import ms_sim_sites

# class Hierarchy(object):
#
#     def __init__(self,H1,H2):
#         pass

def factorial(k,start = 1):
    #compute k * k-1 * ... * start
    assert isinstance(k,int)
    assert k >= 0
    return reduce(lambda x,y: x*y, xrange(start,k+1),1)

def binom(n,k):
    #compute binomial coefficients
    a = factorial(n,n-k+1)
    b = factorial(k)
    #assert a%b == 0
    return a/b

def odd_factorial(n):
    '''returns 1 * 3 * ... * 2n -1'''
    assert isinstance(n,int)
    assert n > 0
    return reduce(lambda x,y: x*y, xrange(1,2*n,2), 1 )

def hist_single_root_with_n_leaves(n):
    #return binom(n,2) * binom(n-1,2) * ... * binom(2,2)
    return reduce( lambda x,y: x*binom(y,2), xrange(2,n+1), 1)

def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    ## copied from https://docs.python.org/2/library/itertools.html
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

def unordered_rooted_leaf_labelled_tree_from_nested_tuples(T):
    if type(T) is tuple:
        subtree_counter = Counter( (unordered_rooted_leaf_labelled_tree_from_nested_tuples(t) for t in T) )
        return LabelledUnorderedRootedTree(subtree_counts = subtree_counter)
    else:
        return LabelledUnorderedRootedTree(rootLabel = T)

def unordered_rooted_leaf_labelled_tree_from_haplotype_matrix(S,seq_labels = None,column_labels = None):

    n,sites = S.shape

    if seq_labels == None:
        seq_labels = np.arange(n)
        #seq_labels = [None]*n

    if column_labels == None:
        column_labels = [None]*sites
        #column_labels = np.arange(sites)

    assert sequences_consistent(S)

    #S1 = sort_matrix_coumns_in_rev_lex_order(S)
    # non_leaf_subtrees = Counter()
    # for i in range(sites):
    #     I1 = np.arange(n)[S1[:,i] == 1] # row indices with mutations in ith collumn
    #     I2 = np.arange(n)[S1[:,i] == 0]# row indices without mutations in ith collumn
    #
    #     T_mut = unordered_rooted_leaf_labelled_tree_from_haplotype_matrix(S1[I1],seq_labels[I1])
    #     non_leaf_subtrees[T_mut] += 1

    # order the columns in lexicographical order
    column_permutation = np.lexsort(np.rot90(S.transpose()))
    #S1 = S[column_permutation]

    #a dictionary mapping each row to the tree associated with that row.
    label2tree = dict( ( (i, LabelledUnorderedRootedTree(rootLabel = seq_labels[i]) ) for i in range(n) ) )
    #tree2label = dict( ( (LabelledUnorderedRootedTree(rootLabel = i), i ) for i in seq_labels ) )
    for c in column_permutation:
        affected_sequences = np.arange(n)[S[:,c] == 1]
        assert len(affected_sequences) < n
        affected_trees = set( (T for s,T in label2tree.iteritems() if s in affected_sequences) )
        T_new = LabelledUnorderedRootedTree(rootLabel = column_labels[c], subtree_counts = Counter(affected_trees))
        for s in affected_sequences:
            label2tree[s] = T_new

    #add the MRCA-root
    trees = set(label2tree.values())
    assert len(trees) > 1
    T = LabelledUnorderedRootedTree(subtree_counts = Counter(trees))

    return T

def sequences_consistent(S):
    # Returns True if all collumn-pairs pass the three gammete test and False
    # otherwise. The ancestral type is presumed to be 0 in every collumn.
    number_of_non_0_elements = np.apply_along_axis(func1d = sum, axis = 0, arr = (S != 0))
    affectedSites = [i for i,n in enumerate(number_of_non_0_elements) if n > 1]
    # discriminant = lambda i: sum(S[j,i] != 0 for j in xrange(S.shape[0])) > 1
    # affectedSites = filter(discriminant, xrange(S.shape[1]))
    for i,s1 in enumerate(affectedSites):
        #for s2 in filter(lambda x: x > s1 , affectedSites):
        for s2 in affectedSites[(i+1):]:
            if not three_gammete_test(S[:,s1], S[:,s2]):
                return False
    return True

def three_gammete_test(c1,c2):
    '''
    Takes two columns/or characters.
    Returns False if the characters contain all three of the gammetes (not 0, 0)
    (not 0, not 0) and (0, not 0). Returns True otherwise.
    '''
    n = len(c1)
    assert n == len(c2)

    AX_occurs = reduce(lambda x,y: x or y, map(lambda i: c1[i] == 0 and c2[i] != 0, range(n)))
    YA_occurs = reduce(lambda x,y: x or y, map(lambda i: c1[i] != 0 and c2[i] == 0, range(n)))
    YX_occurs = reduce(lambda x,y: x or y, map(lambda i: c1[i] != 0 and c2[i] != 0, range(n)))

    return not(AX_occurs and YA_occurs and YX_occurs)

def sort_matrix_coumns_in_rev_lex_order(S):
    return S[:,np.lexsort(np.rot90(S.transpose()))[::-1]]

def labelledUnorderedRootedTree_from_subtrees(subtrees):
    return LabelledUnorderedRootedTree(subtree_counter = Counter(subtrees))

class LabelledUnorderedRootedTree(object):
    '''Represent a rooted unirdered tree ass a collection of subtrees under a labelled root'''

    def __init__(self, rootLabel = None, subtree_counts = Counter() ):
        '''Initialize a tree'''

        #We want trees to be immutable (using them as keys in a hash-table would
        #be ill advised otherwise).
        object.__setattr__(self, "subtree_counts", subtree_counts)
        object.__setattr__(self, "rootLabel", rootLabel)
        object.__setattr__(self, "nodes", self.count_nodes_recursive() )
        object.__setattr__(self, "leaves", self.count_non_root_leaves_recursive() )
        object.__setattr__(self, "rootDegree", self.compute_rootDegree() )
        object.__setattr__(self, "hash_precomputed", self.precompute_hash())

        # self.subtree_counts = subtree_counts
        # self.rootLabel = rootLabel
        # self.nodes = self.count_nodes_recursive()
        # self.leaves = self.count_non_root_leaves_recursive()
        # self.rootDegree = self.compute_rootDegree()

    #To guard against accidentially changing the contents of a tree, we rneder
    #them immutable, by manually making it raise an error to try to modify an
    #instance of this class.
    def __setattr__(self, *args):
        raise TypeError
    def __delattr__(self, *args):
        raise TypeError

    def __hash__(self):
        return self.hash_precomputed

    def precompute_hash(self):
        const = 468395662504823 ## a prime taken from https://primes.utm.edu/curios/page.php/468395662504823.html

        #Hash the root-label
        if self.rootLabel == None:
            rootHash = 0
        else:
            rootHash = hash(self.rootLabel)

        #hash sub-trees
        if sum(self.subtree_counts.values()) == 0:
            return rootHash
        else:

            #For each subtree, we hash that subtree, add one, multiply by our
            #constant, and bit-shift the hashed value by the multiplicity of
            #that subtree. All of the results are then combined as a bitwise
            #xor.
            #Note that since the bitwise xor is commutative, the ordering
            #of subtrees does not matter. This is a feature, not a bug!

            subtree_hash = reduce(lambda x,y: x^y, map(lambda x: ((1+hash(x[0]))*const)<<x[1], self.subtree_counts.items()), 0)

            return rootHash + subtree_hash

    def __eq__(self,other):

        if id(self) == id(other):
            return True

        elif type(self) != type(other):
            return False

        elif self.rootLabel != other.rootLabel:
            return False

        elif hash(self) != hash(other):
            return False

        else:
            return self.subtree_counts == other.subtree_counts

    def count_non_root_leaves_recursive(self):
        ''' returns the number of leaves. if we interpret a rooted tree as having arcs
from a root to the roots of a sub-tree, we want to count the number of leaves
which have out-degree 0 (i.e. the root only counts as a leaf when it has no subtrees)'''
        if len(self.subtree_counts) == 0:
            return 1
        else:
            return sum( [count*T.leaves for T,count in self.subtree_counts.items()] )

    def count_nodes_recursive(self):
        if len(self.subtree_counts) == 0:
            return 1
        else:
            return 1 + sum( [count*T.nodes for T,count in self.subtree_counts.items()] )

    def compute_rootDegree(self):
        return sum(self.subtree_counts.values())

    def __str__(self):
        if self.rootDegree == 0:
            str_root = str(self.rootLabel) if self.rootLabel != None else ''
            return str_root
        else:
            str_root = '%s'%str(self.rootLabel) if self.rootLabel != None else ''
            return '(%s)%s'%(', '.join(map(lambda x: ', '.join([str(x[0])]*x[1]),self.subtree_counts.items())),str_root)

    def getSubtrees(self):
        'returns all subtrees. each subtree has a multiplicity according to its count'
        return sum([ (tree,)*count for tree,count in self.subtree_counts.items()], () )

    def getSubtrees_unique(self):
        return tuple(self.subtree_counts.keys())


def count_histories(tree,computed_histories):

    if tree in computed_histories:
        #Step 1: see if we already know the answer.
        #If we do, there is no reason to redo our work.
        return computed_histories[tree], computed_histories

    if tree.rootDegree == 0:
        #Case: the tree under consideration is a leaf (has no subtrees)
        computed_histories[tree] = 1
        return computed_histories[tree],computed_histories

    if tree.rootDegree == 1:
        #Case: the root has only a single subtree.
        subtree = tree.getSubtrees()[0]
        history_count_new, computed_histories = count_histories(subtree,computed_histories)
        computed_histories[tree] = history_count_new
        return computed_histories[tree], computed_histories

    if max(t.rootDegree for t in tree.getSubtrees_unique() ) == 0:
        #print 'womp!'
        # A special case: T has no propper subtrees, but only leaves below. Here
        # a closed form formula for the number of histories exists, and is much
        # quicker to apply than iterating over all hierarchies.
        computed_histories[tree] = hist_single_root_with_n_leaves(tree.rootDegree)
        return computed_histories[tree], computed_histories

    else:
        subtrees = tree.getSubtrees()
        subtrees_unique = tree.getSubtrees_unique()

        #compute the number of histories of each subtree
        subtree_histories = dict()
        for subtree in subtrees_unique:
            subtree_histories[subtree],computed_histories = count_histories(subtree,computed_histories)

        #print 'subtree_histories = %s'%str(subtree_histories.values())

        #sum over all binary hierarchies over the subtrees
        computed_histories[tree] = reduce(lambda x,H: x + evaluate_hierarchy(H,subtree_histories)[0],  binaryHierarchies(subtrees), 0)

        return computed_histories[tree], computed_histories

def binaryHierarchies(S):
    'returns all labelled binary hierarchies over the set S'

    if len(S) == 1:
        return [( S[0] )]

    #Iterate over all partitions of S into two non-empty subsets S1 and S2. We
    # designate a single element in S, and let it always be a member of S1
    S_list = list(S)
    I = set(range(1,len(S)))
    Hs = []
    for I1 in powerset(range(1,len(S))):
        I1 = set(I1)
        I2 = I.difference(I1)
        assert I1.union(I2) == I
        I1.add(0)
        assert I2.isdisjoint(I1)

        if len(I2) > 0:
            S1 = [S_list[i] for i in I1]
            S2 = [S_list[i] for i in I2]

            #print 'I = %s = I1 + I2 = %s + %s'%(str(I),str(I1),str(I2))

            H1s = binaryHierarchies(S1)
            H2s = binaryHierarchies(S2)

            Hs.extend( ( (H1,H2) for H1 in H1s for H2 in H2s ) )
    return Hs


def reduce_over_tree(H, f_int = lambda x, y: x + y, f_leaf = lambda x: x, f_int_initial = 0):
    ''' take a tree encoded as a tuple of tuples, apply f_leaf to every leaf, and
reduce over f_int at every internal node
e.g.
    reduce_over_tree( ( (a,b), (c,d) ), g, f, f0 )
will return
    f( f( g(a) , f( g(b), f0) ), f( g(c), f(g(d),f0) ) ) '''
    if type(H) == tuple:
        subtree_results = [reduce_over_tree(h, f_int, f_leaf, f_int_initial) for h in H]
        agregated_result = reduce(f_int, subtree_results, f_int_initial)
        #print 'H:\t%s\nsep:\t%s\ncombined:\t%s'%(str(map(str,H)),str(subtree_results),str(agregated_result))
        return agregated_result
        #print map(str,H),subtree_results
    else:
        return f_leaf(H)

def evaluate_hierarchy(Hierarchy,subtree_history_counts):

    def f_leaf(T):
        return (subtree_history_counts[T], T.nodes)

    def f_int(x,y):
        H1,H2 = x[0],y[0]
        K1,K2 = x[1],y[1]
        assert K1>=0
        assert K2>=0
        binomial_foctor = binom(max(K1 + K2 - 2,0), max(K1 - 1,0))
        return (H1 * H2 * binomial_foctor, K1 + K2)

    return reduce_over_tree(Hierarchy, f_int, f_leaf, f_int_initial = (1,0))

# def flatten(H):
#     '''Flatten a hierarchy.
#     i.e. flatten( ( (1,2), (3,4) ) ) = (1,2,3,4)'''
#     if type(H) == tuple:
#         return sum( ( flatten(h) for h in H ), () )
#     else:
#         return (H,)

# def evaluate_hierarchy(Hierarchy,subtree_history_counts):
#     if len(Hierarchy) == 2:
#
#         H1 = evaluate_hierarchy(Hierarchy[0],subtree_history_counts)
#         H2 = evaluate_hierarchy(Hierarchy[1],subtree_history_counts)
#
#         S1 = flatten(Hierarchy[0])
#         S2 = flatten(Hierarchy[1])
#
#         K1 = sum(T.nodes for T in S1) # should not be leaves!
#         K2 = sum(T.nodes for T in S2)
#
#         binomial_foctor = binom(K1 + K2 - 2, K1 - 1)


def testcase_hein_et_al():
    '''Verifies that the tree associated with the dataset
seq1: 1 1 0 0
seq2: 1 1 0 1
seq3: 0 0 1 0
seq4: 0 0 1 0
seq5: 0 0 0 0

has exactly 71 ancestral histories when sites are unlabelled and sequences are labelled'''

    #The dataset corresponds to the tree { { { 1 , {2} } } , {3, 4} , {5} }
    T = unordered_rooted_leaf_labelled_tree_from_nested_tuples( (((1,(2,)),),(3,4),5) )

    print 'considering tree T = %s\n...\n'%str(T),
    hist, hist_dict = count_histories(T,dict())

    try:
        assert hist == 71
        print 'Passed!'
    except AssertionError:
        print 'Filed: incorrect number of histories'
    except:
        print 'Failed: histories could not be computed!'

def testcase_single_root(n_max = 10,supress_individual_tests = False):
    print 'Considering trees  T = {1,2}, {1,2,3}, ... , {1,...,%i}'%n_max
    try:
        for n in range(2,n_max+1):
            T = unordered_rooted_leaf_labelled_tree_from_nested_tuples(tuple(range(n)))
            hist, hist_dict = count_histories(T,dict())
            if supress_individual_tests:
                print '.',
            else:
                print '\n{1,...,%i}\t:\t%s'%(n,str(hist)),
            assert hist == hist_single_root_with_n_leaves(n)
            if not supress_individual_tests: print '  (pass)',
        print '\nAll cases passed!'
    except AssertionError:
        print 'Filed: incorrect number of histories'
    except:
        print 'Failed: histories could not be computed!'


def simulate_data_and_count_histories(n,s,N,seedval = None, print_results = True, add_header = True):

    if seedval != None:
        np.random.seed(seedval)
    #simulate N datasets with n sequences and s seggregating sites.

    simulation = ms_sim_sites(n,s,N,with_trees = False)
    S_list = simulation['S_list']

    #the ms-call that produced the experiment. Kept for reproducibility-purposes.
    command = simulation['metadata']['command']

    #compute the associated trees
    T_list = [unordered_rooted_leaf_labelled_tree_from_haplotype_matrix(S) for S in S_list]

    #compute the number of histories of each dataset
    known_hist = dict()
    hist_list = [count_histories(T,known_hist)[0] for T in T_list]

    #Format each individual tree into a line.
    #(format is csv with semicolons as separators)
    out_strings = ['%i ;\t%i ;\t%i ;\t%s ;\t%s'%(n, s, hist_list[i], str(T_list[i]), command) for i in range(N)]

    if add_header:
        header = 'sequences ; seg sites ;  histories; gene_tree ; command'
        out_strings = [header] + out_strings

    results = '\n'.join(out_strings)

    if print_results:
        print results

    return results

# def hist_recursion(trees,trees_subset):
#     assert trees_subset.issubset(trees)
#     events_total = sum(T.leaves for T in trees)
#     events_subset = sum(T.leaves for T in trees_subset)
#     H_1 =
#     H_2 =
#     coefficient =
#     return
