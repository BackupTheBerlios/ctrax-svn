import kcluster
import kcluster0
import numpy as num
import numpy.random as random
import hungarian

DTYPE = num.float32
state = random.get_state()

def initialize():
    N = 1000
    d = 2
    k = 3
    x = random.randint(low=0,high=100,size=(N,d))
    w = random.uniform(0,100,N)
    return (x.astype(DTYPE),k,w.astype(DTYPE))

def test():
    state = random.get_state()
    (x,k,w) = initialize()
    (mu,S,priors,gamma,err) = kcluster.gmm(x,k,weights=w,nreplicates=10)
    #print "mu = " + str(mu)
    #print "idx = " + str(idx)

#def test1():
#    (x,c) = initialize()
#    D = kcluster.clusterdistfun1(x,c)

def test0():
    random.set_state(state)
    (x,k,w) = initialize()
    (mu,S,priors,gamma,err) = kcluster0.gmm(x,k,weights=w,nreplicates=10)
    #print "mu = " + str(mu)
    #print "idx = " + str(idx)

if __name__ == "__main__":
    for trial in range(10):
        print "trial = " + str(trial)
        (x,k,w) = initialize()
        state = random.get_state()
        if True:
            (mu0,S0,priors0,gamma0,err0) = kcluster0.gmm(x,k,weights=w,nreplicates=10)
            random.set_state(state)
            (mu,S,priors,gamma,err) = kcluster.gmm(x,k,weights=w,nreplicates=10)            
            D = (num.tile(mu[:,0].reshape((k,1)),(1,k)) - mu0[:,0])**2 + \
                (num.tile(mu[:,1].reshape((k,1)),(1,k)) - mu0[:,1])**2
            (a1,a2) = hungarian.hungarian(D)
            if True or \
                    num.max(num.abs(mu - mu0[a1,:])) > .001 or \
                    num.max(num.abs(S-S0[:,:,a1])) > .001 or \
                    num.max(num.abs(priors - priors0[a1])) > .001 or \
                    num.max(num.abs(gamma-gamma0[:,a1])) > .001 or \
                    num.abs(err - err0) > .001:
                print "cluster approx cluster0["+str(a1)+"]"
                print "max(|mu - mu0|) = " + str(num.max(num.abs(mu - mu0[a1,:])))
                print "max(|S - S0|) = " + str(num.max(num.abs(S-S0[:,:,a1])))
                print "max(|priors - priors0|) = " + str(num.max(num.abs(priors - priors0[a1])))
                print "max(|gamma - gamma0|) = " + str(num.max(num.abs(gamma-gamma0)))
                print "(err - err0)/err0 = " + str((err-err0)/err0)

        else:

            (mu,S,priors) = kcluster.gmminit(x,k,weights=w)
            random.set_state(state)
            (mu0,S0,priors0) = kcluster0.gmminit(x,k,weights=w)
            if True or \
                    num.max(num.abs(mu - mu0)) > .001 or \
                    num.max(num.abs(S-S0)) > .001 or \
                    num.max(num.abs(priors - priors0)) > .001:
                print "max(|mu - mu0|) = " + str(num.max(num.abs(mu - mu0)))
                print "max(|S - S0|) = " + str(num.max(num.abs(S-S0)))
                print "max(|priors - priors0|) = " + str(num.max(num.abs(priors - priors0)))
            for i in range(100):
                [gamma,newe] = kcluster.gmmmemberships(mu,S,priors,x,w)
                [gamma0,newe0] = kcluster0.gmmmemberships(mu0,S0,priors0,x,w)
                if True or \
                        num.max(num.abs(gamma-gamma0)) > .001 or \
                        num.abs(newe-newe0) > .001:
                    print "i = %d: max(|gamma - gamma0|) = "%i + str(num.max(num.abs(gamma-gamma0)))
                    print "|newe - newe0| = " + str(num.abs(newe-newe0))
                kcluster.gmmupdate(mu,S,priors,gamma,x,w)
                kcluster0.gmmupdate(mu0,S0,priors0,gamma0,x,w)
                if True or \
                        num.max(num.abs(mu - mu0)) > .001 or \
                        num.max(num.abs(S-S0)) > .001 or \
                        num.max(num.abs(priors - priors0)) > .001:
                    print "i = %d: max(|mu - mu0|) = "%i + str(num.max(num.abs(mu - mu0)))
                    print "max(|S - S0|) = " + str(num.max(num.abs(S-S0)))
                    print "max(|priors - priors0|) = " + str(num.max(num.abs(priors - priors0)))


