import numpy as num
from numpy.random import rand
import numpy.linalg
import scipy.linalg.decomp as decomp

n = 5
d = 2
nclusts = 3

x = rand(n,d)
c = rand(nclusts,d)
#S = rand(d,d,nclusts)
#S = S + S.swapaxes(0,1)
S = num.zeros((d,d,nclusts))
for j in range(nclusts):
    tmp = rand(d,d)
    S[:,:,j] = num.dot(tmp,tmp.T)

print "x = " + str(x)
print "c = " + str(c)
print "S = ["
for j in range(nclusts):
    print " " + str(S[:,:,j])
print "]"

normal = (2.0*num.pi)**(num.double(d)/2.0)

D = num.zeros((nclusts,n))
E = num.zeros((nclusts,n))
for i in range(nclusts):
    D[i,:] = num.sum( (x - c[i,:])**2, axis=1 )
    E[i,:] = num.sum((x - num.tile(c[i,:],[n,1]))**2,axis=1)
    print "D[%d] = "%i + str(D[i,:])
    print "E[%d] = "%i + str(E[i,:])

gamma1 = num.zeros((n,nclusts))
gamma2 = num.zeros((n,nclusts))
for j in range(nclusts):
    print "j = " + str(j)
    print "c.shape = " + str(c.shape)
    diffs = x - c[j,:]
    zz = S[0,0,j]*S[1,1,j] - S[0,1,j]**2
    temp1 = (diffs[:,0]**2*S[1,1,j]
            - 2*diffs[:,0]*diffs[:,1]*S[0,1,j]
            + diffs[:,1]**2*S[0,0,j]) / zz

    print "temp1 = " + str(temp1)
    ch = decomp.cholesky(S[:,:,j])
    temp2 = num.transpose(num.linalg.solve(num.transpose(ch),num.transpose(diffs)))
    temp2 = num.sum(temp2**2,axis=1)
    gamma1[:,j] = num.exp(-.5*temp1)/(normal*num.sqrt(zz))
    gamma2[:,j] = num.exp(-.5*temp2)/(normal*num.prod(num.diag(ch)))
    print "temp2 = " + str(temp2)
    print "sigma1 = " + str(num.sqrt(zz))
    print "sigma2 = " + str(num.prod(num.diag(ch)))
    print "gamma1 = " + str(gamma1[:,j])
    print "gamma2 = " + str(gamma2[:,j])
