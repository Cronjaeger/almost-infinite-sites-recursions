'''
Implementing some recursions for determining the state space size of labelled
versions of the infinite sites model.

Mathias C. Cronjager June 2018
'''

def fallingFactorial(n,k):
    '''compute n x (n-1) x ... x (n-k + 1)'''
    return reduce(lambda x,y: x*y, (n - i for i in xrange(k)), 1)

def factorial(n):
    return fallingFactorial(n,n-1)

def binomial(n,k):
    return fallingFactorial(n,k)/factorial(k)

def multinomial(n, nlist):
    assert n == sum(nlist)
    return reduce(lambda x,y: x/y, map(factorial,nlist) ,factorial(n))

def divisors(n):
    return (x for x in xrange(1,n+1) if n%x == 0)


def memoize(f):
    '''Memoize a function
    Returns a version of f(...) which stores all computed values for later
    reference. The first time f(x) is called, we evaulate f(x) as usual. We then
    proceed to store the key-value-pair (x,f(x)) in a dictionary D. The next
    time f(x) is called, we notice that 'x' is in the keys of D and we just
    retrieve D[x] and return it rather than calling f(x) again.'''

    memory = {}

    def memoized_function(*args):
        try:
            return memory[args]
        except KeyError:
            f_x = f(*args)
            memory[args] = f_x
            return memory[args]

    return memoized_function

@memoize
def t_p_LU(n,l):
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
        b = t_np_LU(n-1,l-1)
        result = a + b
        return result

@memoize
def t_np_LU(n,l):
    if l == n-1:
        return 1
    if 2 <= l < n-1 and n > 2:
        summand = 0
        for n1 in xrange(2,n):
            for l1 in xrange(1,min(l,n1)):
                n2 = n - n1
                l2 = l - l1
                labelFactor = binomial(l,l1)
                t1Count = t_p_LU(n1,l1 + 1) + t_np_LU(n1,l1)
                t2Count = t_p_LU(n2,l2 + int(n2 > 1)) + t_np_LU(n2,l2)
                summand += labelFactor * n2 * t1Count * t2Count
        try:
            assert summand % (n - 1) == 0
        except AssertionError:
            print 'Count (= %i) is not divisible by n-1 (= %i)'%(summand,n-1)
            print ' Make sure the numebrs and boundary-conditions are right!'
        return summand / (n - 1)
    else:
        return 0

def generateStateSpaceTable(nMax = 20, sMax = 5, t_function = t_np_LU):
    header = '\\begin{tabular}{%s}'%('r'*(sMax + 2))
    header += '\n'
    header += 'Sample Size & \\multicolumn{%i}{c}{Segregating Sites}'%(sMax + 1)
    header += '\\\\\n'
    header += ' &' + ' &'.join([str(s) for s in range(0,sMax+1)])
    lines = [header]
    for n in range(2,nMax+1):
        line = '%i '%n
        for s in range(0,sMax+1):
            N = n + s + 1
            l = n
            states = t_function(N,l)
            line += '& {:,}'.format(states)
            # print '%i Lab. seq., %i Ulab. pos : %i'%(n,s,t_np_LU(N,l))
        lines.append(line)
    return '\\\\\n'.join(lines) + '\n\\end{tabular}'

if __name__ == '__main__':
    print generateStateSpaceTable()
