import numpy as num
from ellipses import *
from pylab import *
from kcluster import gmm
import scipy.ndimage.measurements as meas

class ShapeParams:
    def __init__(self,major=0,minor=0,area=0,ecc=0):
        self.major = major
        self.minor = minor
        self.area = area
        self.ecc = ecc

class EllipseParams:
    def __init__(self):
        self.maxshape = ShapeParams(55,25,55*25*num.pi*4)
        self.minshape = ShapeParams(15,5,15*5*num.pi*4)
        self.minbackthresh = .1
        self.maxpenaltymerge = 40
        self.maxareadelete = 5

def ell2cov(a,b,theta):
    S = num.zeros((2,2))
    costheta = num.cos(theta)
    sintheta = num.sin(theta)
    a = a**2
    b = b**2
    S[0,0] = costheta**2*a + sintheta**2*b
    S[1,1] = sintheta**2*a + costheta**2*b
    S[0,1] = sintheta*costheta*(a-b)
    S[1,0] = S[0,1]
    return S

def cov2ell(S):
    # compute axis lengths from covariance matrix
    tmp1 = S[0,0] + S[1,1]
    tmp2 = num.sqrt(4.0*S[0,1]**2.0 + (S[0,0] - S[1,1])**2.0)
    eigA = (tmp1+tmp2)/2.0
    eigB = (tmp1-tmp2)/2.0
    # compute angle
    angle = 0.5*num.arctan2( 2.0*S[0,1], S[0,0] - S[1,1] )
    if eigB > eigA:
        sizeW = num.sqrt(eigA)
        sizeH = num.sqrt(eigB)
    else:
        sizeW = num.sqrt(eigB)
        sizeH = num.sqrt(eigA)
    return (sizeH,sizeW,angle)

def drawellipse(ellipse,format='w',params={}):
    theta = num.linspace(-.03,2*num.pi,100)
    x = 2*ellipse.major*num.cos(theta)
    y = 2*ellipse.minor*num.sin(theta)
    X = num.cos(ellipse.angle)*x - num.sin(ellipse.angle)*y
    Y = num.sin(ellipse.angle)*x + num.cos(ellipse.angle)*y
    X += ellipse.center.x
    Y += ellipse.center.y
    h = plot(X,Y,format,**params)
    return h

def ellipsepixels(ellipse,bounds):
    # convert axes to covariance matrix
    S = ell2cov(ellipse.major,ellipse.minor,ellipse.angle)
    # get all pixels in box
    [x,y] = num.meshgrid(num.arange(bounds[2],bounds[3],1),num.arange(bounds[0],bounds[1],1))
    # subtract off center
    x -= ellipse.center.x
    y -= ellipse.center.y
    # compute Mah distance
    Sinv = num.linalg.inv(S)
    d = x**2*Sinv[0,0] + 2*Sinv[0,1]*x*y + y**2*Sinv[1,1]
    # threshold at 4
    bw = d <= 4
    return bw

def copyellipse(ellipses,i,newellipse):
    ellipses[i].center.x = newellipse.center.x
    ellipses[i].center.y = newellipse.center.y
    ellipses[i].major = newellipse.major
    ellipses[i].minor = newellipse.minor
    ellipses[i].angle = newellipse.angle
    ellipses[i].area = newellipse.area


def weightedregionpropsi(BWI,w):
    # normalize weights
    w = w / sum(w)
    # compute mean
    [r,c] = num.where(BWI)
    centerX = sum(c*w)
    centerY = sum(r*w)
    # compute variance
    r = r - centerY
    c = c - centerX
    S = num.zeros((2,2))
    S[0,0] = sum(w*c**2)
    S[1,1] = sum(w*r**2)
    S[0,1] = sum(w*c*r)
    S[1,0] = S[0,1]
    (sizeH,sizeW,angle) = cov2ell(S)
    # if there is only one pixel in this connected component,
    # then the variance will be 0
    if sizeH < .125 or num.isnan(sizeH):
        sizeH = .125
        sizeW = .125
    elif sizeW < .125 or num.isnan(sizeW):
        sizeW = .125
    area = num.pi * sizeW * sizeH * 4
    return Ellipse( centerX, centerY, sizeW, sizeH, angle, area )

def weightedregionprops(L,ncc,dfore):
    # allocate
    ellipses = EllipseList()
    S = num.zeros((2,2))
    
    # loop through the connected components
    for i in range(ncc):
        BW = L == i+1
        ellipse = weightedregionpropsi(BW,dfore[BW])
        ellipses.append(ellipse)
    return ellipses

def getboundingboxbig(ellipse,sz,params):
    # returns a box that has width and height 2*max major axis length
    # around the center of the ellipse
    r1 = floor(ellipse.center.y-params.maxshape.major*4)
    if r1 < 0: r1 = 0
    r2 = ceil(ellipse.center.y+params.maxshape.major*4)+1
    if r2 > sz[0]: r2 = sz[0]
    c1 = floor(ellipse.center.x-params.maxshape.major*4)
    if c1 < 0: c1 = 0
    c2 = ceil(ellipse.center.x+params.maxshape.major*4)+1
    if c2 > sz[1]: c2 = sz[1]
    return (r1,r2,c1,c2)

def getboundingboxtight(ellipse,sz):
    # returns the bounding box of the ellipse
    r1 = floor(ellipse.center.y-ellipse.major*2)
    if r1 < 0: r1 = 0
    r2 = ceil(ellipse.center.y+ellipse.major*2)+1
    if r2 > sz[0]: r2 = sz[0]
    c1 = floor(ellipse.center.x-ellipse.major*2)
    if c1 < 0: c1 = 0
    c2 = ceil(ellipse.center.x+ellipse.major*2)+1
    if c2 > sz[1]: c2 = sz[1]
    return (r1,r2,c1,c2)

def getnewlabel(Lnewbox,ncc,Lbox,i):
    # choose the component corresponding to this one
    if ncc == 1:
        # if there is only one connected component, then no need to
        # do anything complicated
        llowerthresh = 1
    else:
        # if pixels are connected at a higher threshold, they are also
        # connected at a lower threshold. so choose the label of any
        # member of the original connected component
        newl = Lnewbox[Lbox==i+1]
        llowerthresh = newl[0]
        if num.all(newl==llowerthresh) == False:
            print "Something is wrong -- this should never happen!\n"
            bins = linspace(-.5,ncc+.5,ncc+2)
            votes = num.histogram(Lnewbox[Lbox==i+1],bins)
            llowerthresh = argmax(votes)
    return llowerthresh

def trylowerthresh(ellipses,i,L,dfore,params):
    # try lowering the threshold around the small ellipses
    (r1,r2,c1,c2) = getboundingboxbig(ellipses[i],L.shape,params)
    dforebox = dfore[r1:r2,c1:c2]
    isforebox = dforebox >= params.minbackthresh
    Lbox = L[r1:r2,c1:c2]
    # label connected components in this new foreground image
    (Lnewbox,ncc) = meas.label(isforebox)
    inew = getnewlabel(Lnewbox,ncc,Lbox,i)
    # see if this caused the current connected component to be merged with
    # other connected components
    tmp = Lbox[Lnewbox==inew]
    if num.any(num.logical_and(tmp != 0,tmp != i+1)):
        # if there is a merge, then try merging instead
        return (True,ellipses[i])
    # get properties of this connected component
    ellipsenew = weightedregionpropsi(Lnewbox==inew,dforebox[Lnewbox==inew])
    # see if it is now big enough
    issmall = ellipsenew.area < params.minshape.area
    if issmall == False:
        print 'new area of ellipse is %f > %f' % (ellipsenew.area,params.minshape.area)
        # success!
        copyellipse(ellipses,i,ellipsenew)
    return (issmall,ellipsenew)

def findclosecenters(ellipses,i,params):
    # maximum distance between centers
    maxdmergecenter = params.maxshape.major*4+ellipses[i].minor*2
    
    # indices other than i
    isotherind = num.ones(len(ellipses),dtype=bool)
    isotherind[i] = False
    for j in range(len(ellipses)):
        if ellipses[j].area == 0:
            isotherind[j] = False
    otherinds = num.where(isotherind)[0]
    nother = len(otherinds)

    # first threshold x distance
    dx = num.zeros(nother)
    for j in range(nother):
        dx[j] = abs(ellipses[otherinds[j]].center.x-ellipses[i].center.x)
    isclose = dx < maxdmergecenter
    cannotmerge = num.any(isclose) == False

    # then threshold y distance
    if cannotmerge == False:
        dy = num.zeros(nother)
        dy[:] = maxdmergecenter
        dy[isclose] = abs(ellipses[j].center.y-ellipses[i].center.y)
        isclose[isclose] = dy[isclose] < maxdmergecenter
        cannotmerge = num.any(isclose) == False

    # then threshold Euclidean distance
    if cannotmerge == False:
        maxdmergecentersquared = maxdmergecenter**2
        d = num.zeros(nother)
        d[:] = maxdmergecentersquared
        d[isclose] = dx[isclose]**2 + dy[isclose]**2
        isclose[isclose] = d[isclose] < maxdmergecentersquared
        cannotmerge = num.any(isclose) == False

    indsmerge = otherinds[isclose]
    return indsmerge

def computemergepenalty(ellipses,i,j,L,dfore,params):
    # compute parameters of merged component
    BWmerge = num.logical_or(L == i+1,L == j+1)
    ellipsemerge = weightedregionpropsi(BWmerge,dfore[BWmerge])
    # see if the major, minor, area are small enough
    if (ellipsemerge.area > params.maxshape.area) or (ellipsemerge.minor > params.maxshape.minor) or (ellipsemerge.major > params.maxshape.major):
        return (params.maxpenaltymerge,ellipses[i])
    # find pixels that should be foreground according to the ellipse parameters
    (r1,r2,c1,c2) = getboundingboxtight(ellipsemerge,L.shape)
    isforepredmerge = ellipsepixels(ellipsemerge,num.array([r1,r2,c1,c2]))
    # pixels that were foreground
    isforepredi = ellipsepixels(ellipses[i],num.array([r1,r2,c1,c2]))
    isforepredj = ellipsepixels(ellipses[j],num.array([r1,r2,c1,c2]))
    isforepredi = num.logical_or(isforepredi, (L[r1:r2,c1:c2]==i+1))
    # pixels that are now foreground that weren't before
    newforemerge = num.logical_and(isforepredmerge,num.logical_or(isforepredi,isforepredj)==False)
    # compute the total background score for this new region that must be foreground
    dforemerge = dfore[r1:r2,c1:c2]
    dforemerge = 1 - dforemerge[newforemerge]
    dforemerge[dforemerge<0] = 0
    mergepenalty = num.sum(dforemerge)
    return (mergepenalty,ellipsemerge)

def mergeellipses(ellipses,i,j,ellipsemerge,issmall,L,params):
    # set i to the new merged ellipse
    copyellipse(ellipses,i,ellipsemerge)
    # remove j
    ellipses[j].area = 0
    issmall[i] = ellipsemerge.area < params.minshape.area
    issmall[j] = False
    L[L==j+1] = i+1

def trymerge(ellipses,issmall,i,L,dfore,params):
    # find connected components whose centers are at most maxdmergecenter
    # from the target
    closeinds = findclosecenters(ellipses,i,params)
    
    # if there are no close centers, just return
    if len(closeinds) == 0:
        return False
    # compute the penalty for each close center
    mergepenalty = ones(len(closeinds))
    mergepenalty[:] = params.maxpenaltymerge
    ellipsesmerge = EllipseList()
    for j in range(len(closeinds)):
        (mergepenalty[j],newellipse) = computemergepenalty(ellipses,i,closeinds[j],L,dfore,params)
        print "penalty for merging with ellipse:"
        printellipse(ellipses[closeinds[j]])
        print '= %f' % mergepenalty[j]
                     
        ellipsesmerge.append(newellipse)

    bestjmerge = num.argmin(mergepenalty)
    minmergepenalty = mergepenalty[bestjmerge]

    # see if the penalty is small enough, if not, return
    if minmergepenalty > params.maxpenaltymerge:
        print "merge penalty too large"
        return False

    # perform the merge
    canmergewith = closeinds[bestjmerge]
    mergeellipses(ellipses,i,canmergewith,ellipsesmerge[bestjmerge],issmall,L,params)
    print 'merged %d with %d' % (i,canmergewith)
    return True

def trydelete(ellipses,i,issmall,params):
    if ellipses[i].area < params.maxareadelete:
        ellipses[i].area = 0
        issmall[i] = False

def deleteellipses(ellipses):
    i = 0
    while True:
        if i >= len(ellipses):
            break
        if ellipses[i].area == 0:
            ellipses.pop(i)
        else:
            i+=1

def printellipse(ellipse):
    print '[x: %f y: %f a: %f b: %f t: %f A: %f]' % (ellipse.center.x,ellipse.center.y,ellipse.major,ellipse.minor,ellipse.angle,ellipse.area)
    
def fixsmall(ellipses,L,dfore,params):
    issmall = num.zeros(len(ellipses),dtype=bool)
    
    for i in range(len(ellipses)):
        issmall[i] = ellipses[i].area < params.minshape.area
    while num.any(issmall):
        i = num.where(issmall)[0]
        i = i[0]
        didmerge = False
        print "trying to fix ellipse %d: " % i
        printellipse(ellipses[i])
        (issmall[i],ellipselowerthresh) = trylowerthresh(ellipses,i,L,dfore,params)
        if issmall[i] == False:
            print "Succeeded by lowering threshold:"
            printellipse(ellipselowerthresh)
            
        if issmall[i]:
            print "Could not lower threshold. Trying to merge"
            didmerge = trymerge(ellipses,issmall,i,L,dfore,params)
            print "After attempting to merge, ellipse is now:"
            printellipse(ellipses[i])
            print "didmerge = "
            print didmerge
        if issmall[i] and (didmerge==False):
            print "Could not merge. Trying to delete."
            # set ellipses[i] to be ellipselowerthresh
            copyellipse(ellipses,i,ellipselowerthresh)
            trydelete(ellipses,i,issmall,params)
            print "After deleting, ellipse is:"
            printellipse(ellipses[i])
            issmall[i] = False
    deleteellipses(ellipses)

def normpdfln2x2cov(x,y,S):
    # exponent term
    invS = num.linalg.inv(S)
    p = - (x**2*invS[0,0] + y**2*invS[1,1] + 2*x*y*invS[0,1]) * .5
    # normalization term
    p -= log( 2.0*num.pi*num.sqrt(S[0,0]*S[1,1] - S[0,1]**2) )
    return p

def trysplit(ellipses,i,isdone,L,dfore,params):

    # get datapoints
    (r,c) = num.where(L==i+1)
    x = num.hstack((c.reshape(c.size,1),r.reshape(r.size,1)))
    w = dfore[L==i+1]
    w = w / num.mean(w)
    ndata = r.size
    print 'ndata = %d' % ndata

    # compute the BIC for one cluster
    S0 = ell2cov(ellipses[i].major,ellipses[i].minor,ellipses[i].angle)
    mu0 = num.array([ellipses[i].center.x,ellipses[i].center.y])
    print 'loglik: '
    print normpdfln2x2cov(c-ellipses[i].center.x,r-ellipses[i].center.y,S0)
    loglik = num.sum(w*normpdfln2x2cov(c-ellipses[i].center.x,r-ellipses[i].center.y,S0))
    BIC0 = -2*loglik + num.log(ndata)*5*1
    print 'BIC0 = %f' % BIC0

    # try splitting into more clusters
    ncomponents = 2
    while True:
        (mu,S,priors,gamma,negloglik) = gmm(x,ncomponents,weights=w)
        print 'for ncomponents = %d:' % ncomponents
        print 'mu = '
        print mu
        print 'S = '
        for j in range(ncomponents):
            print S[:,:,j]
        print 'priors = '
        print priors
        print 'negloglik = %f' % negloglik
        BIC = 2*negloglik + num.log(ndata)*5*ncomponents
        if BIC >= BIC0:
            break
        ncomponents += 1
        mu0 = mu
        S0 = S
        gamma0 = gamma
        BIC0 = BIC

    ncomponents -= 1
    if ncomponents == 1:
        isdone[i] = True
    else:
        # get id
        idx = num.argmax(gamma0,axis=1)
        # replace
        ellipses[i].center.x = mu0[0,0]
        ellipses[i].center.y = mu0[0,1]
        (ellipses[i].major,ellipses[i].minor,ellipses[i].angle) = cov2ell(S0[:,:,0])
        # add new
        for j in range(1,ncomponents):
            (a,b,angle) = cov2ell(S0[:,:,j])
            area = a*b*angle*4
            ellipse = Ellipse( mu[j,0], mu[j,1], b, a, angle, area )
            ellipses.append(ellipse)
            L[r[idx==j],c[idx==j]] = len(ellipses)
        isdone.append(num.zeros(ncomponents-1,dtype=bool))
        
def fixlarge(ellipses,L,dfore,params):

    # whether or not we have tried to fix the ellipse
    isdone = num.zeros(len(ellipses),dtype=bool)
    # set to True if the ellipse is not large
    for i in range(len(ellipses)):
        isdone[i] = ellipses[i].area <= params.maxshape.area
    while True:
        # find an ellipse that is not done
        i = num.where(isdone)[0]
        # if there aren't any, break
        if i.size == 0:
            break
        i = i[0]
        trysplit(ellipses,i,isdone,L,dfore,params)

# create parameters
params = EllipseParams()

# make two ellipses
ellipse1 = Ellipse(100,100,20,50,num.pi/3.0)
ellipse2 = Ellipse(200,120,20,50,num.pi/6.0)

# get an image of the two ellipses
bounds = num.array([0,200,0,300])
BW1 = ellipsepixels(ellipse1,bounds)
BW2 = ellipsepixels(ellipse2,bounds)
BW = num.logical_or(BW1,BW2)

# create a noisy version for experimenting with
dfore = num.double(BW)

## experiment with lowering threshold:
#dfore[BW1] = .3
#(r,c) = num.where(BW1)
#dfore[r[0],c[0]] = 1

# experiment with merging
#dfore += dfore*(num.random.standard_normal(size=BW.shape)/2)
#dfore[dfore<0] = 0
#
#BWapprox = dfore >= .5
#
# find connected components
#(Lapprox,nccapprox) = meas.label(BWapprox)
#
# fit ellipses
#ellipsesfit = weightedregionprops(Lapprox,nccapprox,dfore)
#
# display ellipses
#print 'ncc = %d\n' % nccapprox
#for i in range(nccapprox):
#    print 'ellipse %d: [%f, %f, %f, %f, %f].\n' % (i,ellipsesfit[i].center.x,ellipsesfit[i].center.y,ellipsesfit[i].major,ellipsesfit[i].minor,ellipsesfit[i].angle)

## try lowering the threshold
#(issmall,ellipselowerthresh) = trylowerthresh(ellipsesfit,0,Lapprox,dfore,params)
#print 'ellipselowerthresh: [%f, %f, %f, %f, %f].\n' % (ellipselowerthresh.center.x,ellipselowerthresh.center.y,ellipselowerthresh.major,ellipselowerthresh.minor,ellipselowerthresh.angle)

# try merging

#fixsmall(ellipsesfit,Lapprox,dfore,params)
#
#print "after fixing small: "
#print 'ncc = %d\n' % len(ellipsesfit)
#for i in range(len(ellipsesfit)):
#    printellipse(ellipsesfit[i])

# try splitting
L = num.uint16(BW)
ncc = num.max(L)
ellipses = weightedregionprops(L,ncc,dfore)
isdone = num.zeros(len(ellipses),dtype=bool)
trysplit(ellipses,0,isdone,L,dfore,params)

im = imshow(BW, cmap=cm.gray)
for i in range(len(ellipses)):
    h = drawellipse(ellipses[i],format='w',params={"linewidth": 2})
show()
