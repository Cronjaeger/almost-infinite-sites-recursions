from math import log10,log
import numpy as np
from formula_outline import prob_External

def harmonic(n):
    '''Returns the n-th harmonic number.'''
    if n<1:
        return 0
    else:
        return sum([1.0/i for i in xrange(1,n+1)])

def waterson(k,n):
    '''Returns watersons estimator of the muation-rate.
    k -- number of segregating sites
    n -- sample-size
    '''
    if n == 0:
        return 0
    else:
        return k/harmonic(n)

def maxWithIndex(l):
    '''Takes a list l, and returns (max(l), min{i : l(i) == max(l)})'''
    iMax = 0
    maxSoFar = l[0]
    for i,x in enumerate(l):
        if x > maxSoFar:
            maxSoFar = x
            iMax = i
    return maxSoFar,iMax

def logMaximizer(objective_function,
                x0,
                gridpoints = 10,
                steps = 10,
                log10Diameter0 = 2,
                scaleFactor = 2.0/9.0,
                verbose = False
                ):
    """
    Takes some objective-function f: [0,\infty] --> \mathbb{R}, and attempts
    to determine argmax(f), by maximizing over a sequence of logarithmically
    spaced grids, centered about the argMax-value from the last step.
    """
    log_x0 = log10(x0)
    log10Diameter = log10Diameter0
    log_x_max = log_x0
    x_max = x0

    for i in xrange(steps):

        if verbose:
            print "Performing step %i:"%i
            print " x_max = %s"%str(x_max)
            print " log_x_max = %s"%str(log_x_max)
            print " log10Diameter = %s"%str(log10Diameter)

        logGridMin = log_x_max - (log10Diameter / 2.0)
        logStepSize = log10Diameter / (gridpoints - 1)

        if verbose:
            print " logGridMin = %s"%str(logGridMin)
            print " logStepSize = %s"%str(logStepSize)

        logGrid = [logGridMin + i*logStepSize for i in xrange(gridpoints + 1)]
        grid = [10**log_x for log_x in logGrid]

        if verbose:
            print " logGrid:\n  %s"%str(logGrid)
            print " Grid:\n  %s"%str(grid)

        f_values = [ objective_function(x) for x in grid]

        if verbose:
            print " (x,f(x))-pairs :\n  %s"%str(zip(grid,f_values))

        fx_max,i_max = maxWithIndex(f_values)
        x_max,log_x_max = grid[i_max],logGrid[i_max]

        if verbose:
            print " i_max = %i"%i_max
            print " x_max = %s"%str(x_max)
            print " fx_max = %s"%str(fx_max)
            print " log_x_max = %s"%str(log_x_max)
            print "\n"

        log10Diameter *= scaleFactor

    if verbose:
        print " final step comleted."

    return x_max,fx_max

def thetaMLE(
            S,
            nR,
            nC,
            extra_mutations_allowed = 2,
            P = (np.ones((4,4)) - np.eye(4))/3,
            # gridpoints = 6,
            # scaleFactor = 0.5,
            # log10Diameter0 = 2.0,
            # steps = 10,
            verbose = False
            ):

    n = sum(nR)
    deviants = set((1,2,3))
    min_mutations = sum( [ len(deviants.intersection(set([S[i,j] for i in xrange(S.shape[0])]))) * nC[j] for j in xrange(S.shape[1]) ] )
    b = min_mutations + extra_mutations_allowed

    #objective = lambda theta: prob_External(S,nR,nC,b,theta,returnTable = False, P = P)
    def objective(theta):
        return prob_External(S,nR,nC,b,theta,returnTable = False, P = P)

    theta0 = waterson(min_mutations,n)

    if verbose:
        print " n = %i"%n
        print " min_mutations = %i"%min_mutations
        print " b = %i"%b
        print " theta0 = %f"%theta0


    if theta0 == 0.0:
        if verbose:
            print "theta0 == 0 artificially using value 0.001 * ln(n)"
        theta0 = 0.1 * log(n)

    gridpoints = 10
    scaleFactor = 2.0/9.0
    log10Diameter0 = 2.0
    steps = 10

    if verbose:
        print "Running MLE with the following parameters:"
        print " theta0:\t %f"%theta0
        print " steps:\t %i"%steps
        print " gridpoints per step:\t %i"%gridpoints
        print " initial log10Diameter:\t %f"%log10Diameter0
        print " log10Diameter scale-factor per step:\t %f"%scaleFactor

    theta_hat = logMaximizer(objective,theta0,gridpoints,steps,log10Diameter0,scaleFactor,verbose)

    return theta_hat

    #logTheta0 = log10(theta0)
