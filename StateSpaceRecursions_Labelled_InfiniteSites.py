"""
Implementing some recursions for determining the state space size of labelled
versions of the infinite sites model.

Mathias C. Cronjager June 2018
"""

from math import log10, floor
from countingTrees import t_np, t_p, t_s_p, t_s_np

def falling_factorial(n, k):
    # type: (int, int) -> int
    """compute n x (n-1) x ... x (n-k + 1)"""
    return reduce(lambda x, y: x*y, (n - i for i in xrange(k)), 1)


def factorial(n):
    # type: (int) -> int
    return falling_factorial(n, n - 1)


def binomial(n, k):
    # type: (int, int) -> int
    return falling_factorial(n, k) / factorial(k)


def multinomial(n, nlist):
    # type: (int, tuple) -> int
    assert n == sum(nlist)
    return reduce(lambda x,y: x/y, map(factorial, nlist), factorial(n))


def divisors(n):
    return (x for x in xrange(1, n + 1) if n % x == 0)


def memoize(f):
    """Memoize a function
    Returns a version of f(...) which stores all computed values for later
    reference. The first time f(x) is called, we evaulate f(x) as usual. We then
    proceed to store the key-value-pair (x,f(x)) in a dictionary D. The next
    time f(x) is called, we notice that 'x' is in the keys of D and we just
    retrieve D[x] and return it rather than calling f(x) again."""

    memory = {}

    def memoized_function(*args,**kwargs):
        kwargs_tuple = tuple(sorted(kwargs.items()))
        args_all = args + kwargs_tuple
        try:
            return memory[args_all]
        except KeyError:
            f_x = f(*args, **kwargs)
            memory[args_all] = f_x
            return memory[args_all]

    return memoized_function

@memoize
def t_UU(n,l):
    return t_np(long(n) ,long(l), long(1)) if n-l > 1 else 1

@memoize
def t_p_LU(n, l):
    '''count the number of planted unordered trees with n nodes, and l leaves, where
    leaves (nodes with no children) are labelled.'''
    if n==1 and l==1:
        return 1
    if n < 1:
        return 0
    if 2 > l or l > n-1:
        return 0
    else:
        a = t_p_LU(n - 1, l)
        b = t_np_LU(n-1, l-1)
        result = a + b
        return result


@memoize
def t_np_LU(n, l):
    if l == n-1:
        return 1
    if 2 <= l < n-1 and n > 2:
        sum_total = 0
        for n1 in xrange(2, n):
            for l1 in xrange(1, min(l, n1)):
                n2 = n - n1
                l2 = l - l1
                label_factor = binomial(l, l1)
                t1_count = t_p_LU(n1, l1 + 1) + t_np_LU(n1, l1)
                t2_count = t_p_LU(n2, l2 + int(n2 > 1)) + t_np_LU(n2, l2)
                sum_total += label_factor * n2 * t1_count * t2_count
        try:
            assert sum_total % (n - 1) == 0
        except AssertionError:
            print 'Count (= %i) is not divisible by n-1 (= %i)'%(sum_total, n - 1)
            print ' Make sure the numbers and boundary-conditions are right!'
        return sum_total / (n - 1)
    else:
        return 0

@memoize
def t_UL(n,l, k = 0, label_root = False):

    # Verify inputs
    assert n > 0
    assert l > 0
    assert k >= 0
    assert isinstance(label_root, bool)

    if k > 0:
        if n-k == 1 and l == 2: # case: t is a chain, i.e. t = (1)--(2)--...--(n)
            return 1

        elif n-k > 1:
            ways_to_label_chain_below_root = binomial(n - l + int(label_root), k + int(label_root) )
            ways_to_pick_trree_below_chain = t_UL(n-k, l-1, 0)
            return ways_to_label_chain_below_root * ways_to_pick_trree_below_chain
        else:
            return 0

    elif k == 0:
        if n==l==1: # case: t consists of a single root-node
            return 1

        elif n-l == 1 and l>1: # case: t consits of n-1 nodes united directly below a common root.
            return 1

        elif n-l > 1 and l>1: # case: t has at least one internal node, and the root has degree at least two.

            sum_big_t2_terms = 0
            for n1 in xrange(2,n - 1):
                n2 = n - n1
                assert 1 < n2 < n - 1
                assert n1 + n2 == n

                for l1 in xrange(max(l - n2, 1), min(l, n1)):
                    l2 = l-l1
                    assert 1 <= l1 <= n1
                    assert 1 <= l2 <= n2

                    s  = n - l - 1
                    s1 = n1 - l1 - 1 # seg sites on t1 (including root if labelled)
                    s2 = n2 - l2 # seg sites on t2
                    assert s1 >= 0
                    assert s2 >= 0 # if s2==0, t2 has no labels; this case is handled separately (below)
                    assert s == s1 + s2

                    ways_to_pick_root_label = (n - l)**int(label_root)
                    ways_to_partition_non_root_labels = binomial(s, s1)

                    ways_to_pick_tree1 = t_UL(n1, l1, 0, label_root=False)
                    for k1 in xrange(1,n1 - l1 + 2):
                        ways_to_pick_tree1 += t_UL(n1, l1 + 1, k1, label_root=False)

                    ways_to_pick_tree2 = t_UL(n2, l2, 0, label_root=True)
                    for k2 in xrange(1,n2 - l2 + 2):
                        ways_to_pick_tree2 += t_UL(n2, l2 + 1, k2, label_root=True)

                    term  = n2
                    term *= ways_to_pick_root_label
                    term *= ways_to_partition_non_root_labels
                    term *= ways_to_pick_tree1
                    term *= ways_to_pick_tree2

                    sum_big_t2_terms += term
                    assert n1 + n2 == n

            sum_trivial_t2_terms = 0
            for l1 in range(1,l):
                j = l - l1
                n1 = n - j

                ways_to_pick_root_label = (n - l)**int(label_root)

                ways_to_pick_tree1 = t_UL(n1,l1,0)
                for k1 in xrange(1, n1 - l1 + 2):
                    ways_to_pick_tree1 += t_UL(n1, l1 + 1, k1, label_root=False)

                term  = 1
                term *= ways_to_pick_root_label
                term *= ways_to_pick_tree1

                sum_trivial_t2_terms += term

            sum_total = sum_trivial_t2_terms + sum_big_t2_terms
            assert sum_total % (n-1) == 0
            return sum_total // (n-1)

        else:
            return 0
    else:
        raise RuntimeError('k=%i, but should have been handled in previous cases.'%k)

# def t_UL(n , l, k = 0, label_root = False):
#
#     # verify inputs
#     assert n > 0
#     assert l > 0
#     assert k >= 0
#     assert isinstance(label_root, bool)
#
#     if n-k == 1 and (l == 2 or n == l == 1):
#         return 1
#     elif k > 0 and n >= l > 1:
#         if n - k > 1:
#             a = n - l + int(label_root)  # total number of sites to be labelled
#             b = k - 0 + int(label_root)  # total number of sites on chain from root to next node of degree != 2
#             ways_to_pick_labels_to_go_directly_below_root = binomial(a, b)
#             ways_to_pick_tree_below_chain = t_UL(n - k, l - 1, 0, label_root=False)
#             assert ways_to_pick_labels_to_go_directly_below_root >= 0
#             assert ways_to_pick_tree_below_chain >= 0
#             return ways_to_pick_labels_to_go_directly_below_root * ways_to_pick_tree_below_chain
#         else:
#             msg = 'Not sure this should occur: t_UL(%i,%i,%i) was called (returning value 0)' % (n, k, l)
#             raise RuntimeWarning(msg)
#             return 0
#     elif k == 0 and n > l:
#         if n-l == 1 and n>2:
#             return 1
#         sum_total = 0
#         for n2 in xrange(2,n-1):
#             n1 = n - n2
#             assert 2 <= n1 <= n-1
#             for l2 in xrange(max(1, l-n1), min(n2+1, l)):
#                 l1 = l - l2
#                 assert 0 <= l1 < l
#
#                 ways_to_pick_root_label = (n - l) ** int(label_root)
#                 ways_to_partition_non_root_site_labels = binomial(n - l - 1, n1 - l1 - 1)
#
#                 assert 0 < n1 < n
#                 assert 0 < l1 < l
#                 ways_to_label_tree1 = t_UL_aux(n1, l1, label_root=False)
#
#                 assert 0 < n2 < n
#                 assert 0 < l2 < l
#                 ways_to_label_tree2 = t_UL_aux(n2, l2, label_root=True)
#
#                 term = n2
#                 term *= ways_to_pick_root_label
#                 term *= ways_to_partition_non_root_site_labels
#                 term *= ways_to_label_tree1
#                 term *= ways_to_label_tree2
#
#                 sum_total += term
#
#         for j in range(1,l):
#
#             ways_to_pick_root_label = (n - l) ** int(label_root)
#             l1 = l - j
#             n1 = n - j
#
#             term = 0
#             for k1 in xrange(0,n-l1+2):
#                 l_new = l1 + int(k1 > 0 or n1 == 1)
#                 assert l1 > 0
#                 term += t_UL(n1, l1, k1, label_root=label_root)
#
#             term *= ways_to_pick_root_label
#             term *= j
#
#             sum_total += term
#
#         # for l1 in xrange(1,l):
#         #     l2 = l - l1
#         #
#         #     ways_to_pick_root_label = (n - l) ** int(label_root)
#         #     ways_to_pick_tree1 = t_UL_aux(n - l2, l1, label_root=False)
#         #
#         #     term = l2 * ways_to_pick_tree1 * ways_to_pick_root_label
#         #
#         #     sum_total += term
#
#         try:
#             assert sum_total % (n-1) == 0
#         except AssertionError:
#             msg = 'Error in case (n = %i, l = %i, k = %i, label_root = %i): sum_total = %i is not divisible by n-1 = %i'%(n,l,k,int(label_root),sum_total,n-1)
#             raise RuntimeError(msg)
#
#         return sum_total // (n-1)
#
#     else:
#         return 0
#
# def t_UL_aux(n1,l1,label_root=True):
#     assert n1 > 0
#     assert l1 > 0
#     assert isinstance(label_root, bool)
#
#     if l1 == 1:
#         return int(n1 >= 1)
#     else:
#         sum_total = 0
#         for k in xrange(0, n1 - l1 + 1):
#             l_new = l1 + int(k > 0) # l1 is only the number of leaves below the root in the old tree; when k>0, the root is also a leaf.
#             term = t_UL(n1, l_new, k, label_root=label_root)
#             sum_total += term
#         return sum_total

@memoize
def t_LL(n, l, k = 0, label_root = False):

    # verify inputs
    assert n > 0
    assert l > 0
    assert k >= 0
    assert isinstance(label_root, bool)

    if n-k == 1 and (l == 2 or n == l == 1):
        return 1
    elif k > 0:
        if n - k > 1:
            a = n - l + int(label_root) # total number of sites to be labelled
            b = k - 0 + int(label_root) # total number of sites on chain from root to next node of degree != 2
            ways_to_pick_labels_to_go_directly_below_root = binomial(a, b)
            ways_to_pick_tree_below_chain = t_LL(n-k, l-1, 0, label_root = False)
            assert ways_to_pick_labels_to_go_directly_below_root >= 0
            assert ways_to_pick_tree_below_chain >= 0
            return ways_to_pick_labels_to_go_directly_below_root * ways_to_pick_tree_below_chain
        else:
            msg = 'Not sure this should occur: t_LL(%i,%i,%i) was called (returning value 0)'%(n,k,l)
            raise RuntimeWarning(msg)
            return 0
    # elif k == 0 and n - l == 1:
    #     return 1
    elif k == 0:
        sum_total = 0
        for n1 in xrange(2,n):
            n2 = n - n1
            for l1 in xrange(max(l - n2, 1), min(l, n1)):

                l2 = l - l1
                # assert n1 - l1 > 0
                # assert n2 - l2 >= 0

                ways_to_partition_leaf_labels = binomial(l, l1)
                ways_to_pick_root_label = (n - l)**int(label_root)
                ways_to_partition_non_root_site_labels = binomial(n - l - 1, n1 - l1 - 1)

                # assert ways_to_partition_leaf_labels > 0
                # assert ways_to_partition_non_root_site_labels > 0

                ways_to_label_tree1 = t_LL_aux(n1, l1, label_root = False)  # we have already picked the root label.
                ways_to_label_tree2 = t_LL_aux(n2, l2, label_root = True)  # root(t2) is either a leaf or internal.

                # assert ways_to_label_tree1 > 0
                # assert ways_to_label_tree2 > 0

                term = n2
                term *= ways_to_partition_leaf_labels
                term *= ways_to_pick_root_label
                term *= ways_to_partition_non_root_site_labels
                term *= ways_to_label_tree1
                term *= ways_to_label_tree2

                sum_total += term
        try:
            assert sum_total % (n-1) == 0
        except AssertionError:
            msg = 'Error in case (n = %i, l = %i, k = %i, label_root = %i): sum_total = %i is not divisible by n-1 = %i'%(n,l,k,int(label_root),sum_total,n-1)
            raise RuntimeError(msg)
        return sum_total // (n-1)

    else:
        msg = 't_LL(%i,%i,%i) was called. This should not occur (returning value 0).'%(n,k,l)
        raise RuntimeWarning(msg)
        return 0

@memoize
def t_LL_aux(n1, l1, label_root = False):

    assert n1 > 0
    assert l1 > 0
    assert isinstance(label_root, bool)

    if l1 == 1:
        return int(n1 >= 1)
    else:
        sum_total = 0
        for k in xrange(0, n1 - l1 + 1):
            l_new = l1 + int(k > 0) # l1 is only the number of leaves below the root in the old tree; when k>0, the root is also a leaf.
            term = t_LL(n1, l_new, k, label_root = label_root)
            sum_total += term
        return sum_total


# # @memoize
# def t_LL(n, l, k = 0, VERBOSE = False):
#     # type: (int, int, int) -> int
#     """count the number of trees with n nodes and, l leaves and a stalk from the root of length k, subject to the
#     following labelling scheme:
#      - leaves are labelled [x1, ..., xl] (not including the root, unless n=l=1, and k=0),
#      - internal nodes (i.e. nodes of degree > 1) are labelled [y1, ..., ys] (there are n - l - int(k>0) of these)"""
#     if n == l == 1 and k == 0:
#         return 1
#     if k == n - 1 > 0 and l == 2:
#         return 1
#     if k > n - 1 or l < 2 or n < 2 or n < l or k < 0:
#         return 0
#     if k > 0:
#         return binomial(n - l, k) * t_LL(n - k, l-1, 0)
#     if k == 0:
#
#         if l == n-1:
#             return 1
#
#         sum_total = 0
#
#         for n1 in range(2, n):
#
#             n2 = n - n1
#
#             term = 0
#
#             for l1 in xrange(1, l):
#                 l2 = l - l1
#                 site_partitions = binomial(n - l - 1, n1 - l1 - 1)
#                 leaf_partitions = binomial(l, l1)
#                 t1_sum = t_LL(n1, l1, 0) + sum(t_LL(n1, l1 + 1, k1) for k1 in range(1, n - l1 + int(l1 == 1) + 1))
#                 t2_sum = t_star_LL(n2, l2, 0) + sum(t_star_LL(n2, l2 + 1, k2) for k2 in range(1, n - l2 + int(l2 == 1) + 1))
#                 term += site_partitions * leaf_partitions * t1_sum * t2_sum
#                 #if VERBOSE: print '(n = %i, l= %i, k = %i, n1 = %i, n2 = %i, l1 = %i, l2 = %i, t1 = %i, t2 = %i, term = %i)'%(n,l,k,n1,n2,l1,l2,t1_sum,t2_sum,term)
#                 if VERBOSE:
#                     print 'n,l,k = %i, %i, %i'%(n,l,k)
#                     print '  n1, l1, t1 = %i, %i, %i'%(n1, l1, t1_sum)
#                     print '  n2, l2, t2 = %i, %i, %i'%(n2, l2, t2_sum)
#             sum_total += n2 * term
#         # print sum_total
#
#         try:
#             assert sum_total % (n - 1) == 0
#         except AssertionError:
#             msg = 'In case n = %i, l=%i, k = %i: Count (= %i) is not divisible by n-1 (= %i)' % (n, l, k, sum_total, n-1)
#             raise RuntimeError(msg)
#
#         return sum_total / (n - 1)
#
#
# # @memoize
# def t_star_LL(n, l, k):
#     if n == l == 1 and k==0:
#         result = 1
#     elif l == 2 and n > 1 and n - k == 1:
#         result = 1
#     elif k > 0 and n > l > 1 and n - k > 1:
#         result = binomial(n - l, k) * t_star_LL(n - k, l-1, 0)
#     elif k==0 and n > l:
#         result = (n-l) * t_LL(n, l, k)
#     else:
#         result = 0
#     try:
#         assert result >= t_LL(n, l, k)
#     except AssertionError:
#         print 'ERROR: t_star_LL(%i,%i,%i) = %i < t_LL(%i,%i,%i) = %i'%(n, l, k, result, n, l, k, t_LL(n, l, k))
#     return result


def generateCSVTable(n_max = 25,s_max = 25, Verbose = True, log = False):

    regularHeader = 'n, s, seqUlab-siteUlab, seqLab-siteUlab, seqUlab-siteLab, seqLab-siteLab'
    logHeader = 'n, s, log10(seqUlab-siteUlab), log10(seqLab-siteUlab), log10(seqUlab-siteLab), log10(seqLab-siteLab)'
    header = logHeader if log else regularHeader

    if Verbose:
        print header

    lines = [header]
    for n in range(2,n_max+1):
        for s in range(s_max+1):
            N = n+s+1

            if log:
                line = '%i, %i, %.10f, %.10f, %.10f, %.10f'%(n, s, log10(t_UU(N, n)), log10(t_np_LU(N, n)), log10(t_UL(N, n)), log10(t_LL(N, n)))
            else:
                line = '%i, %i, %i, %i, %i, %i'%(n, s, t_UU(N, n), t_np_LU(N, n), t_UL(N, n), t_LL(N, n))

            if Verbose:
                print line

            lines.append(line)

    return '\n'.join(lines)

def identity_function(x):
    """A functtion which does nothing; used as default transformation"""
    return x

def generateLatexStateSpaceTable(n_max = 20, s_max = 5, n_step = 1, s_step = 1, t_function = t_np_LU):
    header = '\\begin{tabular}{l%s}'%('r' * (s_max + 1))
    header += '\n'
    header += '\\hline\n'
    header += 'Sample size & \\multicolumn{%i}{c}{Number of segregating sites}'%(s_max + 1)
    header += '\\\\\n\\hline\n'
    header += '  & ' + '& '.join([str(s) for s in range(0, s_max + 1)])
    lines = [header]
    for n in range(2, n_max + 1, n_step):
        line = '%i '%n
        for s in range(0, s_max + 1, s_step):
            N = n + s + 1
            l = n
            states = t_function(N,l)
            line += '& {:,}'.format(int(floor(states)))
            #line += '& {:.2e}'.format(states)
            #line += '& {:.2g}'.format(states)
            # print '%i Lab. seq., %i Ulab. pos : %i'%(n,s,t_np_LU(N,l))
        lines.append(line)
    return '\\\\\n'.join(lines) + '\\\\\n\\hline\n\\end{tabular}'

def generateLatexStateSpaceTable_multifunction(n_max = 20, s_max = 5, n_step = 1, s_step = 1, transformation = identity_function):
    header = '\\begin{tabular}{lll%s}'%('r' * (s_max + 1))
    header += '\n'
    header += '\\hline\n'
    header += 'Sample size & seq. & pos. & \\multicolumn{%i}{c}{Number of segregating sites}'%(s_max + 1)
    header += '\\\\\n\\hline\n'
    header += '  &  &  & ' + '& '.join([str(s) for s in range(0, s_max + 1)])
    lines = [header]
    for n in range(2, n_max + 1, n_step):
        for f,label in ((t_UU, 'U& U'),
                        (t_np_LU, 'L& U'),
                        (t_UL, 'U& L'),
                        (t_LL, 'L& L')):
            if f == t_UU:
                line = '%i & %s'%(n, label)
            else:
                line = '  & %s'%label

            for s in range(0, s_max + 1, s_step):
                N = n + s + 1
                l = n
                states = f(N, l)
                entry = transformation(states)
                line += '& {:,}'.format(entry)
            lines.append(line)
    return '\\\\\n'.join(lines) + '\\\\\n\\hline\n\\end{tabular}'


def diff_LL_LU(n,l):
    return t_LL(n,l) - t_np_LU(n,l)

def relative_diff(x, y):
    if y!= 0:
        return float(abs(x-y))/y
    elif x!= y:
        return float('inf')
    else:
        return 0


def relative_diff_LL_LU(n,l):
    return relative_diff(t_LL(n, l), t_np_LU(n, l))
    # try:
    #     return ( t_LL(n, l) - t_np_LU(n, l) ) / float(t_LL(n, l))
    # except ValueError:
    #     return float('nan')

def relative_diff_LL_UU(n,l):
    return relative_diff(t_LL(n, l), t_UU(n, l))

def relative_diff_LU_UU(n,l):
    return relative_diff(t_np_LU(n, l), t_UU(n, l))

def relative_diff_UL_UU(n,l):
    return  relative_diff(t_UL(n, l), t_UU(n, l))

def log_relative_diff_LL_LU(n,l):
    try:
        return log10( t_LL(n, l) - t_np_LU(n, l) ) - log10(t_LL(n, l))
    except ValueError:
        return float('nan')

def run_tests_LL():
    assert t_LL(1,1,0) == 1
    assert t_LL(2,2,1) == 1
    assert t_LL(3,2,0) == 1
    assert t_LL(4, 2, 0, label_root=True) == 4
    assert t_LL(5, 3, 0) == 6
    assert t_LL(4,3,1, label_root=True) == 1

    t_630 = t_LL(6,3,0)
    assert t_630 == 30

    for n in xrange(3,10):
        for l in range(1,n):
            try:
                assert diff_LL_LU(n,l) >= 0
            except AssertionError:
                msg  = 'removing internal labels leads to more trees in case n,l = %i,%i'%(n,l)
                msg += '  counts are'
                msg += '   t_LL = %i'%t_LL(n,l)
                msg += '   t_LU = %i' % t_np_LU(n, l)
                raise RuntimeError(msg)

def run_tests_UL():
    assert t_UL(1,1,0) == 1
    assert t_UL(2,2,1) == 1
    assert t_UL(3,2,0) == 1
    assert t_UL(4,3,1) == 1

if __name__ == '__main__':
    #run_tests_LL()
    run_tests_UL()

    # for n in range(2,10):
    #     print 'n = %i'%n
    #     for l in range(1,n + 1):
    #         print '\tt_LL, t_star_LL (%i,%i,0) = %i \t %i'%(n,l,t_LL(n, l , 0), t_LL(n, l , 0, label_root = True))
    #     print ''

    # print t_LL(6,3,0)
    # print t_LL_aux( 5, 2, False)

    # print t_LL(2,1,0)
    # print generateLatexStateSpaceTable(t_function = t_UU)
    # print generateStateSpaceTable(t_function = t_np_LU)
    # print generateStateSpaceTable(t_function = t_LL)
    # print generateStateSpaceTable(t_function = diff_LL_LU)
    # print generateStateSpaceTable(t_function=relative_diff_LL_LU)
    # print t_UL(4,2,0, label_root=True)
    # print t_UL(5,2,0,label_root=True)

    # print generateLatexStateSpaceTable_multifunction(102, 50, n_step=10, s_step=10)
    print generateLatexStateSpaceTable_multifunction(102, 100, n_step=10, s_step=10, transformation=lambda x: int(floor(log10(x))))

    # print generateLatexStateSpaceTable(40, 20, lambda n, l: log10(t_UU(n, l)))
    #print generateLatexStateSpaceTable(40, 20, lambda n, l: log10(t_np_LU(n, l)))
    #print '\n'
    #print generateLatexStateSpaceTable(40, 20, lambda n, l: log10(t_UL(n, l)))
    #print '\n'
    #print generateLatexStateSpaceTable(40, 20, lambda n, l: log10(t_LL(n, l)))
    # #
    # print generateLatexStateSpaceTable(100, 100, relative_diff_LU_UU)
    # print generateLatexStateSpaceTable(100, 100, relative_diff_UL_UU)
    # print generateLatexStateSpaceTable(100, 100, relative_diff_LL_UU)

    #print generateCSVTable(150, 150, False, log=True)
