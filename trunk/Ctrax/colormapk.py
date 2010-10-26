#import matplotlib
import numpy as num
#import matplotlib.colors

def jet(m=64):
    if m==1:
        J = num.array([0,0,0])
        return J
    else:
        m = m - 1

    n = int(num.ceil(m/4.))
    u = num.concatenate((num.linspace(1./n,1.,n),
                         num.ones((n-1)),
                         num.linspace(1.,1./n,n)),0)
    g = num.ceil(n/2.) - int((m%4)==1) + num.r_[0:u.size]
    g = g.astype(int)
    r = g + n
    b = g - n
    g = g[g<m]
    r = r[r<m]
    b = b[b>=0]
    J = num.zeros((m,3));
    J[r,0] = u[:r.size];
    J[g,1] = u[:g.size];
    J[b,2] = u[-b.size:];

    J = num.concatenate((num.zeros((1,3)),J),axis=0)
    return J

def colormap_image(imin,cmap=None,cbounds=None):

    # check to make sure that imin has 2 dimensions or less
    assert(imin.ndim <= 2)

    if cmap is None:
        cmap = jet(64)

    # copy the image
    im = imin.astype(num.double)

    # create the cbounds argument
    if cbounds is not None:
        # check the cbounds input
        assert(len(cbounds)==2)
        assert(cbounds[0]<=cbounds[1])
        # threshold at cbounds
        im = im.clip(cbounds[0],cbounds[1])
        
    neginfidx = num.isneginf(im)
    isneginf = num.any(neginfidx)
    posinfidx = num.isposinf(im)
    isposinf = num.any(posinfidx)
    nanidx = num.isnan(im)
    isnan = num.any(nanidx)
        
    # scale the image to be between 0 and cmap.N-1
    goodidx = [num.logical_not(num.logical_or(num.logical_or(neginfidx,posinfidx),nanidx))]
    minv = im[goodidx].min()
    maxv = im[goodidx].max()
    dv = maxv - minv
    if isposinf:
        maxv = maxv + dv * .025
        im[posinfidx] = maxv
        #print "Some elements of image are +infty"
        
    if isneginf or isnan:
        minv = minv - dv * .025
        im[neginfidx] = minv
        im[nanidx] = minv
        #print "Some elements of image are -infty or nan"

    im = (im - minv) / (maxv-minv)
    im *= (float(cmap.shape[0]-1))

    # round
    im = im.round()
    im = im.astype(int)

    # create rgb image
    r = cmap[im,0]
    g = cmap[im,1]
    b = cmap[im,2]
    newshape = r.shape+(1,)
    rgb = num.concatenate((num.reshape(r,newshape),
                           num.reshape(g,newshape),
                           num.reshape(b,newshape)),2)

    # scale to 0 to 255
    rgb *= 255
    rgb = rgb.astype(num.uint8)

    clim = [minv,maxv]

    return rgb,clim
