"""

ExpBGFGModel

Estimates the probability that a pixel is foreground, background
in any non-aligned video for a given experiment. Background subtraction
based modeling will further refine this model to a single video. The
model implemented here is the following. 

Let

l = pixel label in {fore,back}.
x = pixel location in image.
y = pixel color (or patch) in current image.
d = distance to closest object detection.

Then

P_exp(l | x, y, d) \propto p_exp(y | l, x) p_exp(d | l) P_exp(l | x)

assuming y \perp d | l, x  and  d \perp x | l

We model p_exp(y | l, x) as a Gaussian mu(l,x), sigma(l,x). This is
done by est_fg_appearance_marginal and est_bg_appearance_marginal.

We approximate p_exp(d | l) by histogramming with a fixed bin size. 

We approximate P_exp(l | x) by histogramming with a fixed bin size.

"""
import numpy as num
from params import params
import movies
import annfiles as annot
import estconncomps as est
import scipy.ndimage.morphology as morph
import scipy.ndimage.filters as filters

params.interactive = False
DEBUG = True

if DEBUG:
    import pdb
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm


class ExpBGFGModel:

    def __init__( self, movienames, annnames ):
        
        # labeled movies
        # for now, we will assume that every frame of each movie
        # is labeled
        self.movienames = movienames
        self.annnames = annnames

        self.nmovies = len(movienames)
        if len(self.annnames) != self.nmovies:
            print "Number of annotation and movie files must match"
            raise Exception, "Number of annotation and movie files must match"

        # get movie size
        self.movie = movies.Movie( self.movienames[0], False )
        self.nr = self.movie.get_height()
        self.nc = self.movie.get_width()

        # initialize models
        # initialize each model
        self.init_fg_model()
        self.init_bg_model()
        self.init_obj_detection()
        
    def est_marginals(self):

        for i in range(self.nmovies):

            if DEBUG: print "Computing per-movie model for movie %d"%i
            self.compute_models_permovie(i)

        # pool data for foreground model
        if DEBUG: print "Computing foreground appearance marginal"
        self.est_fg_appearance_marginal()

        # pool data for foreground model
        if DEBUG: print "Computing background appearance marginal"
        self.est_bg_appearance_marginal()

        # normalize histograms of distance to object detections
        if DEBUG: print "Computing object detection distance marginal"
        self.est_obj_detection_marginal()

    def est_fg_appearance_marginal(self):

        # radius we end up pooling to
        self.fg_poolradius = num.zeros((self.nr,self.nc))

        pool_data(self.fg_nsamples_px,self.fg_mu_px,
                  self.fg_sigma_px,
                  params.prior_fg_nsamples_pool,
                  params.prior_fg_pool_radius_factor,
                  self.fg_nsamples,self.fg_mu,
                  self.fg_sigma,self.fg_poolradius,
                  self.fg_min_sigma)

        self.fg_log_Z = .5*num.log(2*num.pi) + num.log(self.fg_sigma)

    def est_bg_appearance_marginal(self):

        # radius we end up pooling to
        self.bg_poolradius = num.zeros((self.nr,self.nc))

        pool_data(self.bg_nsamples_px,self.bg_mu_px,
                  self.bg_sigma_px,
                  params.prior_bg_nsamples_pool,
                  params.prior_bg_pool_radius_factor,
                  self.bg_nsamples,self.bg_mu,
                  self.bg_sigma,self.bg_poolradius,
                  self.bg_min_sigma)

        self.bg_log_Z = .5*num.log(2*num.pi) + num.log(self.bg_sigma)    

    def est_obj_detection_marginal(self):

        # normalize histograms
        Z = num.sum(self.obj_detection_dist_counts_fg)
        self.obj_detection_dist_frac_fg = self.obj_detection_dist_counts_fg / Z
        Z = num.sum(self.obj_detection_dist_counts_bg)
        self.obj_detection_dist_frac_bg = self.obj_detection_dist_counts_bg / Z

    def compute_models_permovie( self, i ):
        """
        compute_models_permovie( self, i )
        
        Estimates the models for movie i.
        """

        # open movie
        self.movie = movies.Movie( self.movienames[i], params.interactive )

        # open annotation
        self.trx = annot.AnnotationFile( self.annnames[i], doreadtrx=True )

        # choose frames to learn from: for now, assume all frames are tracked
        self.firstframe = self.trx.firstframetracked
        self.lastframe = self.trx.lastframetracked

        self.framessample = num.unique(num.round(num.linspace(self.firstframe+1,
                                                              self.lastframe,
                                                              params.prior_nframessample)).astype(int))
        if DEBUG: print "Collecting foreground and background data..."

        for j in range(len(self.framessample)):

            # read in the sampled frame
            self.frame = self.framessample[j]
            self.im,self.timestamp = self.movie.get_frame(self.frame)
            self.trxcurr = self.trx.get_frame(self.frame)

            if DEBUG: print "Collecting data for frame %d"%self.frame

            # create an image which is True inside the ellipses, False outside
            self.compute_isfore()

            # update all the models
            self.update_fg_appearance_marginal()
            self.update_bg_appearance_marginal()
            self.update_obj_detection_marginal()
        
    def compute_log_P_fore_given_appearance(self):

        d2 = ((self.im - self.fg_mu) ./ (2*self.fg_sigma))**2
        self.log_P_fore_given_appearance = -d2 - self.fg_log_Z

    def compute_log_P_fore_given_appearance(self):

        d2 = ((self.im - self.fg_mu) ./ (2*self.fg_sigma))**2
        self.log_P_fore_given_appearance = -d2 - self.fg_log_Z

    def compute_isfore(self):

        self.isfore = num.zeros((self.nr,self.nc),dtype=bool)

        # loop through each ellipse
        for ell in self.trxcurr.itervalues():

            # find pixels inside this ellipse
            (r1,r2,c1,c2) = est.getboundingboxtight(ell,(self.nr,self.nc))
            isfore1 = est.ellipsepixels(ell,num.array([r1,r2,c1,c2]))
            self.isfore[r1:r2,c1:c2] = isfore1

        # compute isback by dilating and flipping
        self.isback = morph.binary_dilation(self.isfore,self.dilation_struct) == False
        
    def init_fg_model(self):
        """
        init_fg_model(self)
        Initializes data structures for the foreground model.
        """
        
        self.fg_nsamples_px = num.zeros((self.nr,self.nc))
        self.fg_mu_px = num.zeros((self.nr,self.nc))
        self.fg_sigma_px = num.zeros((self.nr,self.nc))

        self.fg_nsamples = num.zeros((self.nr,self.nc))
        self.fg_mu = num.zeros((self.nr,self.nc))
        self.fg_sigma = num.zeros((self.nr,self.nc))

        # for morphology
        self.dilation_struct = num.ones((2,2),bool)

        self.fg_min_sigma = 1.

    def init_bg_model(self):
        """
        init_bg_model(self)
        Initializes data structures for the background model.
        """
        
        self.bg_nsamples_px = num.zeros((self.nr,self.nc))
        self.bg_mu_px = num.zeros((self.nr,self.nc))
        self.bg_sigma_px = num.zeros((self.nr,self.nc))

        self.bg_nsamples = num.zeros((self.nr,self.nc))
        self.bg_mu = num.zeros((self.nr,self.nc))
        self.bg_sigma = num.zeros((self.nr,self.nc))

        self.bg_min_sigma = 1.

    def update_fg_appearance_marginal(self):
        """
        update_fg_appearance_marginal(self)
        Update the foreground model based on the current frame
        """
        
        # increment number of samples in foreground pixels
        self.fg_nsamples_px[self.isfore] += 1
        # add to mean
        self.fg_mu_px[self.isfore] += self.im[self.isfore]
        # add to sigma
        self.fg_sigma_px[self.isfore] += self.im[self.isfore]**2

    def update_bg_appearance_marginal(self):
        """
        update_bg_appearance_marginal(self)
        Update the background model based on the current frame
        """

        # increment number of samples in background pixels
        self.bg_nsamples_px[self.isback] += 1
        # add to mean
        self.bg_mu_px[self.isback] += self.im[self.isback]
        # add to sigma
        self.bg_sigma_px[self.isback] += self.im[self.isback]**2

    def update_obj_detection_marginal(self):

        # detect objects in the current image
        self.obj_detect()

        # compute distances to detections
        morph.distance_transform_edt(~self.isobj,distances=self.obj_detection_dist)

        # histogram distances for foreground pixels
        counts,edges = num.histogram(self.obj_detection_dist[self.isfore],
                               bins=self.obj_detection_dist_edges)
        self.obj_detection_dist_counts_fg += counts

        # histogram distances for background pixels
        counts,edges = num.histogram(self.obj_detection_dist[self.isback],
                               bins=self.obj_detection_dist_edges)
        self.obj_detection_dist_counts_bg += counts
                      

    def obj_detect(self):

        self.LoGfilter()

    def init_obj_detection(self):

        # parameters
        self.LoG_hsize = 21
        self.LoG_sigma = 5.
        self.LoG_nonmaxsup_sz = 5
        self.LoG_thresh = .5
        self.difftype = 'darkfglightbg'
        self.obj_detection_dist_nbins = 100

        # create LoG filter
        self.make_LoG_filter()

        # output of LoG filtering
        self.LoG_fil_out = num.zeros((self.nr,self.nc))
        # output of nonmaximal suppression
        self.LoG_fil_max = num.zeros((self.nr,self.nc))

        # output of dist transform
        self.obj_detection_dist = num.zeros((self.nr,self.nc))

        # edges used for histogramming distance to object detections. use log-spacing. 
        n1 = int(num.floor(self.obj_detection_dist_nbins/2))
        n2 = self.obj_detection_dist_nbins - n1
        edges1 = num.arange(-.5,n1-.5)
        edges2 = num.unique(num.round(num.exp(num.linspace(num.log(n1+.5),num.log(max(self.nr,self.nc)+.5),n2+1))-.5)+.5)
        self.obj_detection_dist_edges = num.concatenate((edges1,edges2))
        self.obj_detection_dist_nbins = len(self.obj_detection_dist_edges) - 1
        self.obj_detection_dist_centers = \
            (self.obj_detection_dist_edges[1:] + \
                 self.obj_detection_dist_edges[:-1])/2.

        # histograms of distances to object detections for foreground and background
        self.obj_detection_dist_counts_fg = num.zeros(self.obj_detection_dist_nbins)
        self.obj_detection_dist_counts_bg = num.zeros(self.obj_detection_dist_nbins)

    def make_LoG_filter(self):
        
        # Gaussian
        radius = (self.LoG_hsize-1)/2
        [x,y] = num.meshgrid(num.arange(-radius,radius),
                             num.arange(-radius,radius))
        z = -(x**2+y**2) / (2*self.LoG_sigma**2)
        gau = num.exp(z)

        # sum to 1
        Z = num.sum(gau)
        if Z > 0:
            gau = gau / Z

        # Laplacian
        lap = (x**2 + y**2 - 2*self.LoG_sigma**2) / self.LoG_sigma**4

        self.LoG_fil = lap*gau

        # make sum to zero
        self.LoG_fil = self.LoG_fil - num.mean(self.LoG_fil)

        # zero out small values
        absfil = num.abs(self.LoG_fil)
        self.LoG_fil[absfil<.00001*num.max(absfil)] = 0

        # for light flies on a dark background, flip
        if self.difftype == 'lightfgdarkbg':
            self.LoG_fil = -self.LoG_fil
        

    def LoGfilter(self):
        
        # LoG filter
        filters.correlate(self.im.astype(float),
                          self.LoG_fil,
                          output=self.LoG_fil_out,
                          mode='nearest')

        # depending on difftype, we will only look for positive, negative or both values
        if self.difftype == 'other':
            self.LoG_fil_out = num.abs(self.LoG_fil_out)

        # non-maximal suppression + threshold
        filters.maximum_filter(self.LoG_fil_out,size=self.LoG_nonmaxsup_sz,
                               output=self.LoG_fil_max,mode='constant',cval=0)
        self.isobj = num.logical_and(self.LoG_fil_out == self.LoG_fil_max,
                                     self.LoG_fil_out >= self.LoG_thresh)

    
def create_morph_struct(radius):
    if radius <= 0:
        s = False
    elif radius == 1:
        s = num.ones((2,2),bool)
    else:
        s = morph.generate_binary_structure(2,1)
        s = morph.iterate_structure(s,radius)
    return s
    
    
def integral_image_sum(int_a,r0,r1,c0,c1):
    
    s = int_a[r1,c1]
    idx = c0 > 0
    s[idx] -= int_a[r1[idx],c0[idx]-1]
    idx = r0 > 0
    s[idx] -= int_a[r0[idx]-1,c1[idx]]
    idx = num.logical_and(r0 > 0,c0 > 0)
    s[idx] += int_a[r0[idx]-1,c0[idx]-1]

    return s

def pool_data(nsamples_px,mu_px,sigma_px,nsamples_pool,pool_radius_factor,
              nsamples,mu,sigma,poolradius,min_sigma):

    nr = nsamples_px.shape[0]
    nc = nsamples_px.shape[1]

    # precompute integral image
    int_nsamples_px = num.cumsum(num.cumsum(nsamples_px,axis=0),axis=1)
    int_mu_px = num.cumsum(num.cumsum(mu_px,axis=0),axis=1)
    int_sigma_px = num.cumsum(num.cumsum(sigma_px,axis=0),axis=1)

    # whether we are done pooling
    notdone = num.ones((nr,nc),dtype=bool)

    radius = 1
    maxradius = max(nr,nc)
    while (radius <= maxradius) and notdone.any():

        # where do we need to continue incrementing radius?
        r,c = num.where(notdone)

        # current box
        r0 = num.maximum(0,r-radius)
        r1 = num.minimum(nr-1,r+radius)
        c0 = num.maximum(0,c-radius)
        c1 = num.minimum(nc-1,c+radius)

        # compute number of samples in each box
        nsamplescurr = integral_image_sum(int_nsamples_px,r0,r1,c0,c1)

        # which pixels are we now done for?
        newdone = nsamplescurr >= nsamples_pool

        if not newdone.any():
            # increase the radius
            radius = max(radius+1,int(round(radius*pool_radius_factor)))

            if DEBUG: print "no pixels newly with enough samples, incrementing radius to %d"%radius

            continue

        notdone[notdone] = ~newdone

        # for pixels that just finished
        r = r[newdone]
        c = c[newdone]
        r0 = r0[newdone]
        r1 = r1[newdone]
        c0 = c0[newdone]
        c1 = c1[newdone]

        # store the radius and number of samples
        poolradius[r,c] = radius
        nsamples[r,c] = nsamplescurr[newdone]

        # compute the average color
        mu[r,c] = integral_image_sum(int_mu_px,r0,r1,c0,c1) / \
            nsamples[r,c]

        # compute the color standard deviation
        sigma[r,c] = num.sqrt(integral_image_sum(int_sigma_px,r0,r1,c0,c1) / \
            nsamples[r,c])

        if DEBUG: print "n pixels added at radius %d: %d"%(radius,len(r))

        #if DEBUG:
        #    print "\nNew pixels at radius = %d:"%radius
        #    for i in range(len(r)):
        #        print "For pixel (%d,%d), radius = %d, mu = %f, sigma = %f, nsamples = %d"%\
        #          (r[i],c[i],radius,mu[r[i],c[i]],sigma[r[i],c[i]],
        #           nsamples[r[i],c[i]])

        # increase the radius
        radius = max(radius+1,int(round(radius*pool_radius_factor)))

    # make sure standard deviations are large enough
    sigma = num.minimum(sigma,min_sigma)

def test1():

    movienames = ['/media/data/data3/pGAWBCyO-GAL4_TrpA_Rig1Plate01Bowl1_20100820T140408/movie.fmf',]
    annnames = ['/media/data/data3/pGAWBCyO-GAL4_TrpA_Rig1Plate01Bowl1_20100820T140408/movie.fmf.ann',]

    self = ExpBGFGModel(movienames,annnames)
    
    self.est_marginals()

    plt.subplot(221)
    plt.imshow(self.fg_mu,cmap=cm.gray,vmin=0,vmax=255)
    plt.colorbar()
    plt.title("foreground mean")

    plt.subplot(222)
    plt.imshow(self.fg_sigma)
    plt.colorbar()
    plt.title("foreground standard deviation")

    plt.subplot(223)
    plt.imshow(self.bg_mu,cmap=cm.gray,vmin=0,vmax=255)
    plt.colorbar()
    plt.title("background mean")

    plt.subplot(224)
    plt.imshow(self.bg_sigma)
    plt.colorbar()
    plt.title("background standard deviation")

    plt.figure()
    sz = (self.nr,self.nc,1)
    for j in range(4):
        self.frame = self.framessample[j]
        self.im,self.timestamp = self.movie.get_frame(self.frame)
        self.trxcurr = self.trx.get_frame(self.frame)
        self.compute_isfore()
        plt.subplot(1,4,j+1)
        tmp = num.concatenate((self.im.reshape(sz),num.zeros(sz,dtype=num.uint8),255*self.isfore.astype(num.uint8).reshape(sz)),2)
        plt.imshow(tmp)
        plt.title("frame %d"%self.frame)


    plt.figure()
    plt.plot(self.obj_detection_dist_centers,self.obj_detection_dist_frac_bg,'k.-',
             self.obj_detection_dist_centers,self.obj_detection_dist_frac_fg,'r.-')
    plt.legend(('bg','fg'))
    plt.xlabel('Dist to obj detection')
    plt.ylabel('P(dist | label)')

    plt.show()

    return self

def test2():

    movienames = ['/media/data/data3/pGAWBCyO-GAL4_TrpA_Rig1Plate01Bowl1_20100820T140408/movie.fmf',]
    annnames = ['/media/data/data3/pGAWBCyO-GAL4_TrpA_Rig1Plate01Bowl1_20100820T140408/movie.fmf.ann',]

    self = ExpBGFGModel(movienames,annnames)

    # open movie
    i = 0
    self.movie = movies.Movie( self.movienames[i], params.interactive )
    
    # open annotation
    self.trx = annot.AnnotationFile( self.annnames[i], doreadtrx=True )

    self.frame = 1
    self.im,self.timestamp = self.movie.get_frame(self.frame)
    self.trxcurr = self.trx.get_frame(self.frame)
    self.LoGfilter()

    plt.subplot(131)
    plt.imshow(self.LoG_fil_out)
    plt.title('LoG filter')
    plt.colorbar()
    
    plt.subplot(132)
    plt.imshow(self.LoG_fil_out>=self.LoG_thresh)
    plt.title('LoG filter >= thresh')
    plt.colorbar()
    
    plt.subplot(133)
    plt.imshow(self.im,cmap=cm.gray,vmin=0,vmax=255)
    r,c = num.where(self.isobj)
    plt.plot(c,r,'r.')
    plt.axis('image')
    plt.title('isobj')
    plt.show()
    
    return self
    
