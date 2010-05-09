# mt_histgmm.py
# KMB Matlab
# translated to Python 4/12/07 JAB

from matplotlib.mlab import normpdf
import numpy as num
from scipy.linalg.basic import eps, kron # kron substitutes for repmat
import sys

PRINT_PROGRESS = False

def make_vect( arr ):
    return arr.reshape( arr.size )

def histgmm( x, w, maxiters=100, minchange=None, minsigma=None ):
    """Gaussian mixture modeling using 2 clusters.
    Returns 2-D arrays 'mu', 'sigma', 'prior'."""

    if len(w.shape) == 1:
        n = w.size
        m = 1
        w = w.reshape( (n,m) )
        x = x.reshape( (n,m) )
    else:
        n,m = w.shape
    w /= w.sum( 0 )
    x = kron( num.ones((1,m)), x )

    if minchange is None:
        minchange = (x[-1]-x[0]) / n/2
    if minsigma is None:
        minsigma = (x[-1]-x[0]) / n/2

    ## initialize ##

    # first initialization scheme
    mu0 = num.zeros( (2,m) )
    mu0[0,:] = x[0,:]
    mu0[1,:] = x[-1,:]
    sigma0 = (x[-1]-x[0]) / 4 * num.ones( (2,m) )
    gamma = 0.5 * num.ones( (n,m) )
    prior = 0.5 * num.ones( (2,m) )

    mu, sigma, prior, score1 = histgmm_main( mu0, sigma0, gamma, prior, w, x,
                                             maxiters, minchange, minsigma, n, m )

    # second initialization scheme
    c = num.r_[num.zeros(m), w.cumsum()]
    gamma = (c[:-1] < 0.5).astype( num.float )
    gamma.resize( (n,m) )
    mu0 = num.zeros( (2,m) )
    sigma0 = num.zeros( (2,m) )
    tmp = w * gamma
    prior[0,:] = tmp.sum( 0 )
    isprior = prior[0,:] > 0
    # the 'if' statement is necessary because sum( [] ) returns 0, not empty
    if isprior.any():
        mu0[0,isprior] = (tmp * x[:,isprior]).sum(0) / prior[0,isprior]
        a = (x[:,isprior].T - \
              kron( num.ones((n,1)), mu0[0,isprior] )).T
        sigma0[0,isprior] = num.sqrt( (tmp * a**2).sum(0) / prior[0,isprior] )
    tmp = w * (1.-gamma)
    prior[1,:] = tmp.sum( 0 )
    isprior = prior[1,:] > 0
    if isprior.any():
        mu0[1,isprior] = (tmp * x[:,isprior]).sum(0) / prior[1,isprior]
        a = (x[:,isprior].T - \
              kron( num.ones((n,1)), mu0[1,isprior] )).T
        sigma0[1,isprior] = num.sqrt( (tmp * a**2).sum(0) / prior[1,isprior] )
    sigma0[sigma0 < minsigma] = minsigma

    mu2, sigma2, prior2, score2 = histgmm_main( mu0, sigma0, gamma, prior, w, x,
                                                maxiters, minchange, minsigma, n, m )

    ## choose best of mu, mu2
    
    # vectorize scores to allow indexing
    inds2 = num.array((score2,0.)) > num.array((score1,1.))
    if inds2.any():
        mu[:,inds2] = mu2[:,inds2]
        sigma[:,inds2] = sigma2[:,inds2]
        prior[:,inds2] = prior2[:,inds2]

    noprior = prior[0,:] < eps
    if noprior.any():
        mu[0,noprior] = mu[1,noprior]
        sigma[0,noprior] = sigma[1,noprior]
    noprior = prior[1,:] < eps
    if noprior.any():
        mu[1,noprior] = mu[0,noprior]
        sigma[1,noprior] = sigma[0,noprior]

    return mu, sigma, prior

def histgmm_main( mu0, sigma0, gamma, prior, w, x,
                  maxiters, minchange, minsigma, n, m ):
    """Do EM."""
    mu = num.zeros( mu0.shape )
    sigma = num.zeros( mu0.shape )

    ischange = num.ones( m ).astype( num.bool )

    if PRINT_PROGRESS: print "EM iterations:", ; sys.stdout.flush()
    for iter in range( maxiters ):
        if PRINT_PROGRESS: print iter, ; sys.stdout.flush()

        # compute probability for each
        p1 = normpdf( make_vect( x[:,ischange] ).T,
                      make_vect( kron( num.ones((n,1)), mu0[0,ischange] ) ).T,
                      make_vect( kron( num.ones((n,1)), sigma0[0,ischange] ) ) ) \
                      * make_vect( kron( num.ones((n,1)), prior[0,ischange] ) )
        p2 = normpdf( make_vect( x[:,ischange] ).T,
                      make_vect( kron( num.ones((n,1)), mu0[1,ischange] ) ).T,
                      make_vect( kron( num.ones((n,1)), sigma0[1,ischange] ) ) ) \
                      * make_vect( kron( num.ones((n,1)), prior[1,ischange] ) )
        tmp = p1 + p2
        p1[tmp < eps] = 0.5
        p2[tmp < eps] = 0.5
        tmp[tmp < eps] = 1
        gamma[:,ischange] = (p1/tmp).reshape( (n,len(ischange.nonzero()[0])) )

        # update the mean, variance, and prior
        tmp = w[:,ischange] * gamma[:,ischange]
        prior[0,ischange] = tmp.sum( 0 )
        tmp = tmp[:,prior[0,ischange] > 0]
        isprior = prior[0,:] > 0
        c_and_p = num.logical_and( ischange, isprior )
        if c_and_p.any():
            mu[0,c_and_p] = (tmp * x[:,c_and_p]).sum(0) / prior[0,c_and_p]
            a = (x[:,c_and_p].T - \
                  kron( num.ones((n,1)), mu[0,c_and_p] )).T
            sigma[0,c_and_p] = num.sqrt( (tmp * a**2).sum(0) / prior[0,c_and_p] )
        
        tmp = w[:,ischange] * (1-gamma[:,ischange])
        prior[1,ischange] = tmp.sum( 0 )
        tmp = tmp[:,prior[1,ischange] > 0]
        isprior = prior[1,:] > 0
        c_and_p = num.logical_and( ischange, isprior )
        if c_and_p.any():
            mu[1,c_and_p] = (tmp * x[:,c_and_p]).sum(0) / prior[1,c_and_p]
            a = (x[:,c_and_p].T - \
                  kron( num.ones((n,1)), mu[1,c_and_p] )).T
            sigma[1,c_and_p] = num.sqrt( (tmp * a**2).sum(0) / prior[1,c_and_p] )

        sigma[sigma < minsigma] = minsigma

        # see if there is a change
        ischange = num.logical_or( (num.abs( mu - mu0 ) >= minchange),
                                   (num.abs( sigma - sigma0 ) >= minchange) ).any(0)
        if not ischange.any(): break

        mu0 = mu.copy()
        sigma0 = sigma.copy()
    # endfor: EM iterations
    if PRINT_PROGRESS: print

    p1 = normpdf( make_vect( x ).T,
                  make_vect( kron( num.ones((n,1)), mu[0,:] ) ).T,
                  make_vect( kron( num.ones((n,1)), sigma[0,:] ) ) ) \
                  * make_vect( kron( num.ones((n,1)), prior[0,:] ) )
    p2 = normpdf( make_vect( x ).T,
                  make_vect( kron( num.ones((n,1)), mu[1,:] ) ).T,
                  make_vect( kron( num.ones((n,1)), sigma[1,:] ) ) ) \
                  * make_vect( kron( num.ones((n,1)), prior[1,:] ) )
    p1.reshape( (n,m) )
    p2.reshape( (n,m) )
    a = make_vect(w) * (make_vect(gamma)*p1 + (1-make_vect(gamma))*p2)
    score = a.sum( 0 )

    return mu, sigma, prior, score
