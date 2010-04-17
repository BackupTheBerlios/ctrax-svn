import numpy as num
import scipy.cluster.vq as vq
cimport numpy as num
cimport cython

DTYPE = num.float64
ctypedef num.float64_t DTYPE_t
cdef unsigned int d = 2
cdef bint DEBUG = False
cdef bint DEBUG_REPEATABLE_BEHAVIOR = False
#from version import DEBUG, DEBUG_REPEATABLE_BEHAVIOR
cdef DTYPE_t MINPRIOR = .01

cdef inline DTYPE_t square(DTYPE_t x):
    return x*x

@cython.boundscheck(False) # turn of bounds-checking for entire function
def clusterdistfun(num.ndarray[DTYPE_t,ndim=2] x not None,
                   num.ndarray[DTYPE_t,ndim=2] c not None):
    assert x.dtype == DTYPE and c.dtype == DTYPE
    cdef unsigned int nclusts = c.shape[0]
    cdef unsigned int i
    cdef num.ndarray[DTYPE_t,ndim=2] D = num.empty((c.shape[0],x.shape[0]),
                                                   dtype=DTYPE)
    for i in xrange(c.shape[0]):
        D[i,:] = (x[:,0] - c[i,0])**2 + \
            (x[:,1] - c[i,1])**2 # 2d
    return D

@cython.boundscheck(False) # turn of bounds-checking for entire function
def furthestfirst(num.ndarray[DTYPE_t,ndim=2] x not None,
                  unsigned int k, 
                  num.ndarray[DTYPE_t,ndim=1] mu0=None,
                  start='mean'):
    # index
    cdef unsigned int i
    cdef unsigned int j
    # data dimensionality: always 2
    # d = x.shape[1]
    
    # returned centers
    cdef num.ndarray[DTYPE_t,ndim=2] mu = num.empty((k,d),dtype=DTYPE)
    
    # initialize first center
    if mu0 is None:
        if start == 'mean':
            # choose the first mean to be the mean of all the data
            mu[0,:] = num.mean(x,axis=0).astype(DTYPE)
        else:
            # choose a random point to be the first center
            mu[0,:] = x[num.random.randint(0,x.shape[0],1)[0],:]
    else:
        mu[0,:] = mu0
    
    # distance to closest center
    cdef num.ndarray[DTYPE_t,ndim=1] minD = \
        (x[:,0] - mu[0,0])**2 + \
        (x[:,1] - mu[0,1])**2 # 2d
    cdef num.ndarray[unsigned int,ndim=1] idx = num.zeros(x.shape[0],
                                                         dtype=num.uint32)
    cdef num.ndarray[DTYPE_t,ndim=1] newD
    
    for i in xrange(1,k):
        # choose the point furthest from all centers as the next center
        j = num.argmax(minD)
        mu[i,:] = x[j,:]
        # compute distance to current center
        newD = (x[:,0] - mu[i,0])**2 + \
               (x[:,1] - mu[i,1])**2
        # replace distance if this center is closer
        idx[newD < minD] = i
        minD[newD < minD] = newD[newD < minD]

    return (mu,idx)


@cython.boundscheck(False) # turn of bounds-checking for entire function
def gmminit(num.ndarray[DTYPE_t,ndim=2] x not None,
            unsigned int k,
            weights=1,
            unsigned int kmeansiter=20,
            double kmeansthresh=.001):

    cdef num.ndarray[DTYPE_t,ndim=2] mu
    cdef num.ndarray[unsigned int,ndim=1] idx
    cdef num.ndarray[DTYPE_t,ndim=2] D
    cdef num.ndarray[DTYPE_t,ndim=3] S
    cdef num.ndarray[DTYPE_t,ndim=1] priors
    cdef unsigned int i
    cdef num.ndarray[DTYPE_t,ndim=2] diffs
    cdef unsigned int nidx
    cdef bint isweights = isinstance(weights,num.ndarray)
    cdef num.ndarray[DTYPE_t,ndim=2] w

    # initialize using furthest-first clustering
    (mu,idx) = furthestfirst(x,k,start='random')

    # use k-means, beginning from furthest-first
    mu = vq.kmeans(x,mu,kmeansiter,kmeansthresh)[0]

    # get labels for each data point
    D = clusterdistfun(x,mu)
    idx = num.argmin(D,axis=0).astype(num.uint32)

    # allocate covariance and priors
    S = num.empty((d,d,k),dtype=DTYPE)
    priors = num.empty(k,dtype=DTYPE)

    if isweights:
        # unweighted
        for i in xrange(k):
            # compute prior for each cluster
            nidx = num.sum(num.uint32(idx==i))
            priors[i] = nidx
            # compute mean for each cluster
            mu[i,:] = num.mean(x[idx==i,:],axis=0).astype(DTYPE)
            # compute covariance for each cluster
            diffs = x[idx==i,:] - mu[i,:].T
            S[:,:,i] = num.dot(diffs.T,diffs) / nidx
    else:
        w = weights.reshape(x.shape[0],1)
        for i in xrange(k):
            # compute prior for each cluster
            priors[i] = num.sum(w[idx==i])
            # compute mean for each cluster
            mu[i,:] = num.sum(w[idx==i]*x[idx==i,:],axis=0)/priors[i]
            # compute covariance for each cluster
            diffs = x[idx==i,:] - mu[i,:].reshape(1,d)
            diffs *= num.sqrt(w[idx==i])
            S[:,:,i] = num.dot(num.transpose(diffs),diffs) / priors[i]

    # normalize priors
    priors = priors / num.sum(priors)

    # return
    return (mu,S,priors)

@cython.boundscheck(False) # turn of bounds-checking for entire function
def gmm(num.ndarray[DTYPE_t,ndim=2] x,
        unsigned int k,
        weights=None,
        unsigned int nreplicates=10,
        unsigned int kmeansiter=20,
        double kmeansthresh=.001,
        unsigned int emiters=100,
        double emthresh=.001,
        double mincov=.01):

    # for debugging only: reseed the random number generator at 0 for repeatable behavior
    if DEBUG_REPEATABLE_BEHAVIOR:
        num.random.seed(0)

    # number of data points
    cdef unsigned int n = x.shape[0]
    # dimensionality of each data point: 2d
    #d = x.shape[1]
    # initialize min error
    cdef DTYPE_t err
    cdef unsigned int rep
    cdef num.ndarray[DTYPE_t,ndim=2] mubest
    cdef num.ndarray[DTYPE_t,ndim=3] Sbest
    cdef num.ndarray[DTYPE_t,ndim=1] priorsbest
    cdef num.ndarray[DTYPE_t,ndim=2] gammabest
    cdef DTYPE_t minerr = num.inf

    # replicate many times
    for rep in xrange(nreplicates):

        # initialize randomly
        (mu,S,priors) = gmminit(x,k,weights,kmeansiter,kmeansthresh)
        #print "rep = %d"%rep
        #print "mu initialized to " + str(mu)
        #print "S initialized to " + str(S)
        #print "priors initialized to " + str(priors)

        # optimize fit using EM
        (mu,S,priors,gamma,err) = gmmem(x,mu,S,priors,weights,emiters,emthresh,mincov)
        #print "after iterating, "
        #print "1 mu = " + str(mu)
        #print "1 S = " + str(S)
        #print "1 priors = " + str(priors)
        #print "1 gamma = " + str(gamma)
        #print "1 err = " + str(err)

        if rep == 0 or err < minerr:
            mubest = mu
            Sbest = S
            priorsbest = priors
            minerr = err
            gammabest = gamma

    return (mubest,Sbest,priorsbest,gammabest,minerr)

@cython.boundscheck(False) # turn of bounds-checking for entire function
def gmmmemberships(num.ndarray[DTYPE_t,ndim=2] mu,
                   num.ndarray[DTYPE_t,ndim=3] S,
                   num.ndarray[DTYPE_t,ndim=1] priors,
                   num.ndarray[DTYPE_t,ndim=2] x,
                   weights=1,
                   num.ndarray[DTYPE_t,ndim=3] initcovars=None):

    if initcovars is None:
        initcovars = S.copy()

    # number of data points
    cdef unsigned int n = x.shape[0]
    # number of clusters
    cdef unsigned int k = mu.shape[0]

    # allocate output
    cdef num.ndarray[DTYPE_t,ndim=2] gamma = num.empty((n,k),dtype=DTYPE)

    # normalization constant
    cdef double normal = 2*num.pi

    cdef unsigned int j
    cdef num.ndarray[DTYPE_t,ndim=2] diffs
    cdef DTYPE_t zz
    cdef num.ndarray[DTYPE_t,ndim=1] temp

    # compute the activation for each data point
    for j in xrange(k):
        diffs = x - mu[j,:]
        zz = S[0,0,j]*S[1,1,j] - square(S[0,1,j])
        if zz <= 0:
            if DEBUG: print 'S[:,:,%d] = '%j + str(S[:,:,j]) + ' is singular'
            if DEBUG: print 'Reverting to initcovars[:,:,%d] = '%j + str(initcovars[:,:,j])
            S[:,:,j] = initcovars[:,:,j]
            zz = S[0,0,j]*S[1,1,j] - square(S[0,1,j])

        temp = (diffs[:,0]**2*S[1,1,j]
                - 2*diffs[:,0]*diffs[:,1]*S[0,1,j]
                + diffs[:,1]**2*S[0,0,j]) / zz

        gamma[:,j] = num.exp(-.5*temp)/(normal*num.sqrt(zz))

    # include prior
    gamma *= priors
    #print 'after including prior, gamma = '
    #print gamma

    # compute negative log likelihood
    cdef DTYPE_t e
    
    #print "weights = " + str(weights)
    #print "sum(gamma,axis=1) = " + str(num.sum(gamma,axis=1))
    e = -num.sum(num.log(num.sum(gamma,axis=1))*weights)
    #print '1 e = %f' % e

    cdef num.ndarray[DTYPE_t,ndim=1] s
    s = num.sum(gamma,axis=1)
    #print 's = '
    #print s
    # make sure we don't divide by 0
    s[s==0] = 1
    gamma /= s.reshape(n,1)
    #print 'gamma = '
    #print gamma
    
    return (gamma,e)    

@cython.boundscheck(False) # turn of bounds-checking for entire function
def gmmupdate(num.ndarray[DTYPE_t,ndim=2] mu,
              num.ndarray[DTYPE_t,ndim=3] S,
              num.ndarray[DTYPE_t,ndim=1] priors,
              num.ndarray[DTYPE_t,ndim=2] gamma,
              num.ndarray[DTYPE_t,ndim=2] x,
              weights=None,
              double mincov=.01,
              num.ndarray[DTYPE_t,ndim=3] initcovars=None):

    cdef bint isweights = isinstance(weights,num.ndarray)

    if initcovars is None:
        initcovars = S
    
    # number of data points
    cdef unsigned int n = gamma.shape[0]
    # number of clusters
    cdef unsigned int k = gamma.shape[1]
    # dimensionality of data: 2d
    #d = x.shape[1]
    cdef unsigned int i

    if isweights:
        gamma *= weights.reshape(n,1)

    # update the priors (note that it has not been normalized yet)
    priors = num.sum(gamma,axis=0)
    cdef num.ndarray[DTYPE_t,ndim=1] Z = priors.copy()
    cdef DTYPE_t sumpriors = num.sum(priors)
    cdef bint issmall
    if sumpriors > 0:
        priors /= sumpriors
        issmall = num.any(priors < MINPRIOR)
    else:
        if DEBUG: print "All priors are too small, reinitializing"
        if DEBUG: print "gamma = " + str(gamma)
        if DEBUG: print "mu = " + str(mu)
        if DEBUG: print "S = " + str(S)
        issmall = True
        #issmall = num.array((1,)).any() # num.bool_ NOT bool
    #print 'updated priors to ' + str(priors)

    # if any prior is to small, then reinitialize that cluster
    if issmall:
        fixsmallpriors(x,mu,S,priors,initcovars,gamma)
        if DEBUG: print 'directly after fixsmallpriors, priors is ' + str(priors)
        priors = num.sum(gamma,axis=0)
        Z = priors.copy()
        priors /= num.sum(priors)
        priors[priors < MINPRIOR] = MINPRIOR
        priors /= num.sum(priors)
        if DEBUG: print 'after fixsmallpriors, priors is ' + str(priors)

    #issmall = issmall.any()
    if issmall:
        if DEBUG: 
            print 'outside fixsmallpriors'
            print 'reset mu = ' + str(mu)
            for i in xrange(k):
                print 'reset S[:,:,%d] = '%i + str(S[:,:,i])
            print 'reset priors = ' + str(priors)
        #print 'reset gamma = '
        #print gamma

    cdef num.ndarray[DTYPE_t,ndim=2] diffs
    cdef DTYPE_t tr2
    cdef DTYPE_t det
    cdef DTYPE_t mineigval

    for i in xrange(k):
        # update the means
        mu[i,:] = num.sum(gamma[:,i].reshape(n,1)*x,axis=0)/Z[i]
        if DEBUG:
            if issmall: 
                print 'updated mu[%d,:] to '%i + str(mu[i,:])
        # update the covariances
        diffs = x - mu[i,:]
        #if issmall: print 'diffs = ' + str(diffs)
        diffs *= num.sqrt(gamma[:,i].reshape(n,1))
        #if issmall: print 'weighted diffs = ' + str(diffs)
        S[:,:,i] = (num.dot(num.transpose(diffs),diffs)) / Z[i]
        if DEBUG:
            if issmall: 
                print 'updated S[:,:,%d] to [%.4f,%.4f;%.4f,%.4f]'%(i,S[0,0,i],S[0,1,i],S[1,0,i],S[1,1,i])
        # make sure covariance is not too small
        if mincov > 0:
            # hard-coded 2 x 2
            tr2 = (S[0,0,i] + S[1,1,i])/2
            det = S[0,0,i]*S[1,1,i] - S[0,1,i]*S[1,0,i]
            mineigval = tr2 - num.sqrt(square(tr2) - det)

            if num.isnan(mineigval) or mineigval < mincov:
                if DEBUG: 
                    print 'mineigval = %.4f'%mineigval
                    print "tr/2 = " + str(tr2)
                    print "det = " + str(det)
                    print 'S[:,:,%d] = '%i + str(S[:,:,i])
                    print 'reinitializing covariance'
                    print 'initcovars[:,:,%d] = [%.4f,%.4f;%.4f,%.4f]'%(i,initcovars[0,0,i],initcovars[0,1,i],initcovars[1,0,i],initcovars[1,1,i])
                S[:,:,i] = initcovars[:,:,i]

          
@cython.boundscheck(False) # turn of bounds-checking for entire function
def gmmem(num.ndarray[DTYPE_t,ndim=2] x,
          num.ndarray[DTYPE_t,ndim=2] mu0,
          num.ndarray[DTYPE_t,ndim=3] S0,
          num.ndarray[DTYPE_t,ndim=1] priors0,
          weights=1,
          unsigned int niters=100,
          double thresh=.001,
          double mincov=.01):

    #print 'data = '
    #print x
    cdef num.ndarray[DTYPE_t,ndim=2] mu = mu0.copy()
    cdef num.ndarray[DTYPE_t,ndim=3] S = S0.copy()
    cdef num.ndarray[DTYPE_t,ndim=1] priors = priors0.copy()

    cdef DTYPE_t e = num.inf
    cdef DTYPE_t mineigval
    cdef DTYPE_t tr2
    cdef DTYPE_t det
    # store initial covariance in case covariance becomes too small
    if mincov > 0:
        for i in xrange(S.shape[2]):
            #print 'S initially is: '
            #print S[:,:,i]
            if num.any(num.isnan(S[:,:,i])):
                if DEBUG: print "initial S[:,:,%d] contains nans"%i
                if DEBUG: print "S[:,:,%d] = "%i
                if DEBUG: print str(S[:,:,i])
                S[:,:,i] = num.eye(d)*mincov
                
            else:
                tr2 = (S[0,0,i] + S[1,1,i])/2
                det = S[0,0,i]*S[1,1,i] - S[0,1,i]*S[1,0,i]
                mineigval = tr2 - num.sqrt(square(tr2) - det)
                if num.isnan(mineigval) or mineigval < mincov:
                    if DEBUG: print "initial S[:,:,%d] is singular"%i
                    if DEBUG: print "S[:,:,%d] = "%i
                    if DEBUG: print str(S[:,:,i])
                    D, U = num.linalg.eig(S[:,:,i])
                    D[D<mincov] = mincov
                    S[:,:,i] = num.dot(num.dot(U,num.diag(D)),U.transpose())
                    if DEBUG: print 'S[:,:,%d] reinitialized to: '%i
                    if DEBUG: print S[:,:,i]
    cdef num.ndarray[DTYPE_t,ndim=3] initcovars = S.copy()

    cdef unsigned int iter
    cdef num.ndarray[DTYPE_t,ndim=2] gamma
    cdef DTYPE_t newe
    for iter in xrange(niters):

        # E-step: compute memberships
        [gamma,newe] = gmmmemberships(mu,S,priors,x,weights,initcovars)

        # M-step: update parameters
        gmmupdate(mu,S,priors,gamma,x,weights,mincov,initcovars)

        # if we've converged, break
        if newe >= e - thresh:
            break
        
        e = newe

    [gamma,e] = gmmmemberships(mu,S,priors,x,weights,initcovars)

    return (mu,S,priors,gamma,e)
      

@cython.boundscheck(False) # turn of bounds-checking for entire function
def fixsmallpriors(num.ndarray[DTYPE_t,ndim=2] x,
                   num.ndarray[DTYPE_t,ndim=2] mu,
                   num.ndarray[DTYPE_t,ndim=3] S,
                   num.ndarray[DTYPE_t,ndim=1] priors,
                   num.ndarray[DTYPE_t,ndim=3] initcovars,
                   num.ndarray[DTYPE_t,ndim=2] gamma):

    #print 'calling fixsmallpriors with: '
    #print 'mu = ' + str(mu)
    #print 'S = '
    #for i in xrange(S.shape[2]):
    #    print S[:,:,i]
    #print 'priors = ' + str(priors)
    #for i in xrange(initcovars.shape[2]):
    #    print 'initcovars[:,:,%d]: '%i
    #    print initcovars[:,:,i]
    #    print 'initcovars[:,:,%d].shape: '%i + str(initcovars[:,:,i].shape)
        
    cdef num.ndarray[int,ndim=1] smalli
    smalli, = num.where(num.isnan(priors) | (priors < .01))

    if DEBUG:
        print "x = " + str(x)
        print "mu = " + str(mu)
        print "priors = " + str(priors)
        print "smalli = " + str(smalli)

    if smalli.shape[0] == 0:
        return

    # loop through all small priors
    cdef unsigned int i
    cdef num.ndarray[DTYPE_t,ndim=1] p
    cdef unsigned int j
    cdef DTYPE_t newe

    for i in smalli:

        if DEBUG: print 'fixing cluster %d with small prior = %f: '%(i,priors[i])

        # compute mixture density of each data point
        p = num.sum(gamma*priors,axis=1)

        #print 'samples: '
        #print x
        #print 'density of each sample: '
        #print p

        # choose the point with the smallest probability
        j = p.argmin()

        if DEBUG: print 'lowest density sample: x[%d] = '%j + str(x[j,:])

        # create a new cluster
        mu[i,:] = x[j,:]
        S[:,:,i] = initcovars[:,:,i]
        priors *= (1. - MINPRIOR)/(1.-priors[i])
        priors[i] = MINPRIOR

        if DEBUG: 
            print 'reset cluster %d to: '%i
            print 'mu = ' + str(mu[i,:])
            print 'S = '
            print S[:,:,i]
            print 'S.shape: ' + str(S[:,:,i].shape)
            print 'priors = ' + str(priors)

        # update gamma
        [gamma,newe] = gmmmemberships(mu,S,priors,x,1,initcovars)
