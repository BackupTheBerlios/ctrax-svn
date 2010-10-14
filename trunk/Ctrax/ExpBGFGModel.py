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
import bg
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

    def __init__( self, movienames, annnames,
                  LoG_hsize=21, 
                  LoG_sigma=5.,
                  LoG_nonmaxsup_sz=5, 
                  LoG_thresh=.5,
                  difftype='darkfglightbg',
                  obj_detection_dist_nbins=100,
                  isback_dilation_hsize=2,
                  fg_min_sigma = 1.,
                  bg_min_sigma = 1.,
                  min_always_bg_mask_frac = 1.):
        """
        ExpBGFGModel(movienames,annnames)

        Initialize the ExpBGFGModel. 

        Inputs: 
        movienames is a list of all movie files to train on
        annnames is a list of the corresponding annotation files

        The pixel-appearance foreground model is initialized using
        init_fg_model.
        The pixel-appearance background model is initialized using
        init_bg_model.
        The foreground detection model is initialized using
        init_obj_detection.

        Optional inputs:
        LoG_hsize: Width, height of the LoG filter
        (radius = (LoG_hsize-1)/2, so this should be an odd number)
        [21]. 
        LoG_sigma: Standard deviation of the LoG filter [5].
        LoG_nonmaxsup_sz: Size of the neighborhood over which nonmaximal
        suppression will be performed for LoG detection [5]. 
        LoG_thresh: Threshold for LoG object detection [.5].
        difftype: Whether the flies are dark on a light background 
        ('darkfglightbg'), light on a dark background ('lightfgdarkbg'),
        or other ('other') ['darkfglightbg']. 
        obj_detection_dist_nbins: Number of log-spaced bins for distance
        to nearest object detection [100]. 
        isback_dilation_hsize: Size of ones matrix for eroding ~isfore
        to compute isback [2]. 
        fg_min_sigma: Minimum standard deviation for the foreground
        pixel appearance [1]. 
        bg_min_sigma: Minimum standard deviation for the background 
        pixel appearance [1]. 
        
        """

        # parameters
        self.LoG_hsize = LoG_hsize
        self.LoG_sigma = LoG_sigma
        self.LoG_nonmaxsup_sz = LoG_nonmaxsup_sz
        self.LoG_thresh = LoG_thresh
        self.difftype = difftype
        self.obj_detection_dist_nbins = obj_detection_dist_nbins
        self.isback_dilation_hsize = isback_dilation_hsize
        self.fg_min_sigma = fg_min_sigma
        self.bg_min_sigma = bg_min_sigma
        self.min_always_bg_mask_frac = min_always_bg_mask_frac

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
        params.GRID.setsize((self.nr,self.nc))

        # initialize models
        # initialize each model
        self.init_fg_model()
        self.init_bg_model()
        self.init_obj_detection()
        self.init_always_bg_mask()
        
    def est_marginals(self):
        """
        est_marginals()

        This function calls functions to estimate the following 
        marginals from the training data:
        p(pixel appearance | foreground),
        p(pixel appearance | background), 
        p(object detection features | foreground), and
        p(object detection features | background)
        """

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
        
        if DEBUG: print "Estimating always_bg_mask"
        self.est_always_bg_mask()
        
    def est_fg_appearance_marginal(self):
        """
        est_fg_appearance_marginal()

        This function estimates
        p(pixel appearance | foreground)
        from the training data. It uses pool_data to find an estimate
        of the pixel appearance mean and standard deviation for 
        foreground pixels observed near each pixel location. 
        """

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
        """
        est_bg_appearance_marginal()

        This function estimates
        p(pixel appearance | background)
        from the training data. It uses pool_data to find an estimate
        of the pixel appearance mean and standard deviation for 
        background pixels observed near each pixel location. 
        """

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
        """
        est_obj_detection_marginal()

        This function estimates
        p(dist to nearest object detection|foreground) and
        p(dist to nearest object detection|background).

        It uses the histogrammed counts as an estimate. 
        """

        # normalize histograms
        Z = num.sum(self.obj_detection_dist_counts_fg)
        self.obj_detection_dist_frac_fg = self.obj_detection_dist_counts_fg / Z
        Z = num.sum(self.obj_detection_dist_counts_bg)
        self.obj_detection_dist_frac_bg = self.obj_detection_dist_counts_bg / Z

    def compute_models_permovie( self, i ):
        """
        compute_models_permovie( i )
        
        For movienames[i], this function samples frames from throughout
        the movie and computes the data structures necessary for 
        estimating the marginal distributions. 
        """

        # open movie
        self.movie = movies.Movie( self.movienames[i], params.interactive )

        # background model
        self.bg_imgs = bg.BackgroundCalculator(self.movie)

        # open annotation
        self.trx = annot.AnnotationFile( self.annnames[i], bg_imgs=self.bg_imgs, 
                                        doreadtrx=True )

        # choose frames to learn from: for now, assume all frames are tracked
        self.firstframe = self.trx.firstframetracked
        self.lastframe = self.trx.lastframetracked

        self.framessample = num.unique(num.round(num.linspace(self.firstframe+1,
                                                              self.lastframe,
                                                              params.prior_nframessample)).astype(int))
        # update the always_bg mask data structures for this movie
        self.update_always_bg_mask()

        if DEBUG: print "Collecting foreground and background data..."

        for j in range(len(self.framessample)):

            # read in the sampled frame
            self.frame = self.framessample[j]
            self.im,self.timestamp = self.movie.get_frame(self.frame)
            self.trxcurr = self.trx.get_frame(self.frame)

            if DEBUG: print "Collecting data for frame %d"%self.frame

            # create an image which is True inside the ellipses, False outside
            self.isfore_mask()

            # update all the observation models
            self.update_fg_appearance_marginal()
            self.update_bg_appearance_marginal()
            self.update_obj_detection_marginal()
                                
    def compute_log_lik_appearance_given_fore(self):
        """
        compute_log_lik_appearance_given_fore()
        
        computes the log likelihood of each pixel appearance in
        self.im given that the pixel is foreground. 
        """

        d2 = ((self.im - self.fg_mu) / (2*self.fg_sigma))**2
        self.log_lik_appearance_given_fore = -d2 - self.fg_log_Z

    def compute_log_lik_appearance_given_back(self):
        """
        compute_log_lik_appearance_given_back()
        
        computes the log likelihood of each pixel appearance in
        self.im given that the pixel is background. 
        """

        d2 = ((self.im - self.bg_mu) / (2*self.bg_sigma))**2
        self.log_lik_appearance_given_back = -d2 - self.bg_log_Z

    def compute_log_lik_dist_obj(self):
        """
        compute_log_lik_dist_obj()
        
        computes the log likelihood of the distance to the nearest object detection
        for each pixel in the image self.im given that the pixel is in foreground and
        given that the pixel is in background. 
        """
        
        # detect objects in the current image
        self.obj_detect()

        # compute distances to detections
        morph.distance_transform_edt(~self.isobj,distances=self.obj_detection_dist)

        # find which bin each distance falls in
        idx = num.digitize(self.obj_detection_dist,self.obj_detection_dist_edges)

        # lookup probabilities for each bin idx
        self.log_lik_dist_obj_given_fore = num.log(self.obj_detection_dist_frac_fg[idx])
        self.log_lik_dist_obj_given_back = num.log(self.obj_detection_dist_frac_bg[idx])

    def compute_log_lik(self):

        self.compute_log_lik_appearance_given_fore()
        self.compute_log_lik_appearance_given_back()
        self.compute_log_lik_dist_obj()
        self.log_lik_given_fore = self.log_lik_appearance_given_fore + \
            self.log_lik_dist_obj_given_fore
        self.log_lik_given_back = self.log_lik_appearance_given_back + \
            self.log_lik_dist_obj_given_back
            
    def isfore_mask(self):
        """
        isfore_mask()

        Computes a mask self.isfore which is 1 inside of ellipses
        in the current frame and 0 outside. Also computes a mask
        self.isback which is an eroded version of the complement. 
        """

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
        init_fg_model()

        Initializes data structures for the foreground model.
        """
        
        self.fg_nsamples_px = num.zeros((self.nr,self.nc))
        self.fg_mu_px = num.zeros((self.nr,self.nc))
        self.fg_sigma_px = num.zeros((self.nr,self.nc))

        self.fg_nsamples = num.zeros((self.nr,self.nc))
        self.fg_mu = num.zeros((self.nr,self.nc))
        self.fg_sigma = num.zeros((self.nr,self.nc))

        # for morphology
        self.dilation_struct = num.ones((self.isback_dilation_hsize,self.isback_dilation_hsize),bool)


    def init_bg_model(self):
        """
        init_bg_model()

        Initializes data structures for the background model.
        """
        
        self.bg_nsamples_px = num.zeros((self.nr,self.nc))
        self.bg_mu_px = num.zeros((self.nr,self.nc))
        self.bg_sigma_px = num.zeros((self.nr,self.nc))

        self.bg_nsamples = num.zeros((self.nr,self.nc))
        self.bg_mu = num.zeros((self.nr,self.nc))
        self.bg_sigma = num.zeros((self.nr,self.nc))

    def update_fg_appearance_marginal(self):
        """
        update_fg_appearance_marginal()

        Update the foreground pixel appearance data structures 
        based on the current frame. 
        """
        
        # increment number of samples in foreground pixels
        self.fg_nsamples_px[self.isfore] += 1
        # add to mean
        self.fg_mu_px[self.isfore] += self.im[self.isfore]
        # add to sigma
        self.fg_sigma_px[self.isfore] += self.im[self.isfore]**2

    def update_bg_appearance_marginal(self):
        """
        update_bg_appearance_marginal()

        Update the background pixel appearance data structures 
        based on the current frame.
        """

        # increment number of samples in background pixels
        self.bg_nsamples_px[self.isback] += 1
        # add to mean
        self.bg_mu_px[self.isback] += self.im[self.isback]
        # add to sigma
        self.bg_sigma_px[self.isback] += self.im[self.isback]**2

    def update_obj_detection_marginal(self):
        """
        update_obj_detection_marginal()

        Update the object detection data structures based on the
        current frame. 
        """        

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
                      
                      
    def update_priors(self):
        
        self.counts_fg[self.isfore] += 1
        self.counts_tot += 1

    def obj_detect(self):
        """
        obj_detect()

        Detect objects in the current frame. For now, this is just
        LoG filtering. 
        """

        self.LoGfilter()

    def init_obj_detection(self):
        """
        init_obj_detection()

        Initialize object detection data structures. 
        """

        # create LoG filter
        self.make_LoG_filter()

        # output of LoG filtering
        self.LoG_fil_out = num.zeros((self.nr,self.nc))
        # output of nonmaximal suppression
        self.LoG_fil_max = num.zeros((self.nr,self.nc))

        # output of dist transform
        self.obj_detection_dist = num.zeros((self.nr,self.nc))

        # edges used for histogramming distance to object detections. use bins 
        # centered at 0:1:nbins/2, log spacing between nbins/2 and maximum distance 
        # in the image
        n1 = int(num.floor(self.obj_detection_dist_nbins/2))
        n2 = self.obj_detection_dist_nbins - n1
        edges1 = num.arange(-.5,n1-.5)
        edges2 = num.unique(num.round(num.exp(num.linspace(num.log(n1+.5),num.log(num.sqrt(self.nr**2+self.nc**2)+.5),n2+1))-.5)+.5)
        self.obj_detection_dist_edges = num.concatenate((edges1,edges2))
        self.obj_detection_dist_nbins = len(self.obj_detection_dist_edges) - 1
        self.obj_detection_dist_centers = \
            (self.obj_detection_dist_edges[1:] + \
                 self.obj_detection_dist_edges[:-1])/2.

        # histograms of distances to object detections for foreground and background
        self.obj_detection_dist_counts_fg = num.zeros(self.obj_detection_dist_nbins)
        self.obj_detection_dist_counts_bg = num.zeros(self.obj_detection_dist_nbins)

    def init_always_bg_mask(self):
        self.always_bg_mask_count = num.zeros((self.nr,self.nc))
        self.n_always_bg_mask_samples = 0

    def make_LoG_filter(self):
        """
        make_LoG_filter()

        Initialize the LoG filter
        """
        
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
        elif self.difftype == 'darkfglightbg':
            pass
        else:
            self.LoG_fil = num.abs(self.LoG_fil)
        

    def LoGfilter(self):
        """
        LoGfilter()

        Apply LoG filtering to the current frame to detect blobs. 
        """
        
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
                                     
    def update_always_bg_mask(self):

        # TODO: use calibration
        
        # current always bg mask
        always_bg_mask = self.bg_imgs.isarena == False
        
        # count number of times each pixel is always background
        self.always_bg_mask_count += always_bg_mask.astype(float)
        self.n_always_bg_mask_samples += 1

    def est_always_bg_mask(self):
        
        self.always_bg_mask_frac = self.always_bg_mask_count / self.n_always_bg_mask_samples
        self.always_bg_mask = self.always_bg_mask_frac >= self.min_always_bg_mask_frac        
    
def create_morph_struct(radius):
    """
    create_morph_struct(radius)

    Create disk structure for morphological filtering. 
    
    """
    if radius <= 0:
        s = False
    elif radius == 1:
        s = num.ones((2,2),bool)
    else:
        s = morph.generate_binary_structure(2,1)
        s = morph.iterate_structure(s,radius)
    return s
    
    
def integral_image_sum(int_a,r0,r1,c0,c1):
    """
    integral_image_sum(int_a,r0,r1,c0,c1)

    Compute the sum of all values in the box [r0,r1],[c0,c1]
    using the pref-computed integral image int_a. 
    """
    
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

    """
    pool_data(nsamples_px,mu_px,sigma_px,nsamples_pool,pool_radius_factor,
    nsamples,mu,sigma,poolradius,min_sigma)

    At each pixel location, we want to sample nsamples_pool pixels 
    near the current location with the given label. We start by looking
    in a box of radius 1 (so the width and height of the box are 2*1 + 1)
    and count the number of training samples. If this is at least 
    nsamples_pool, then we compute the mean and standard deviation of
    the training samples, and this is the mean and standard deviation
    for the current pixel location. Otherwise, we increase the sample
    radius by the factor pool_radius_factor, and repeat. We continue
    in this way until enough samples are found or until the pooling
    radius exceeds maxradius. 

    Inputs:
    
    nsamples_px is an array of size [nr,nc], and stores the number 
    of training samples of the given label observed at each pixel 
    location. 
    mu_px is an array of size [nr,nc], and stores the mean of all
    the training samples of the given label observed at each pixel 
    location
    sigma_px is an array of size [nr,nc], and stores the std of all
    the training samples of the given label observed at each pixel 
    location
    nsamples_pool is the ideal number of samples we will take around
    each pixel location. 
    pool_radius_factor is the factor to increase the pooling radius
    by at each iteration. 
    min_sigma is the minimum standard deviation to be returned. 

    Modified inputs:

    nsamples is an array of size [nr,nc], and nsamples[i,j] is the
    number of samples collected for pixel location [i,j]. 
    mu is an array of size [nr,nc], and mu[i,j] is the estimated mean
    for pixel location [i,j]
    sigma is an array of size [nr,nc], and sigma[i,j] is the 
    estimated standard deviation for pixel location [i,j]
    poolradius is an array of size [nr,nc], and poolradius[i,j] is
    the radius that samples were taken over for pixel location [i,j]
    """

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
    sigma = num.maximum(sigma,min_sigma)

def test1():
    """
    test1()

    Compute experiment bg and fg models for a set of movies set
    at the start of the function, and plot the results. 
    """

    movienames = ['E:\Data\FlyBowl\Test1\GMR_13B06_AE_01_TrpA_Rig1Plate01BowlA_20101007T161822/movie.ufmf']
    annnames = ['E:\Data\FlyBowl\Test1\GMR_13B06_AE_01_TrpA_Rig1Plate01BowlA_20101007T161822/movie.ufmf.ann',]

    self = ExpBGFGModel(movienames,annnames)
    
    self.est_marginals()

    plt.subplot(231)
    plt.imshow(self.fg_mu,cmap=cm.gray,vmin=0,vmax=255)
    plt.colorbar()
    plt.title("foreground mean")

    plt.subplot(233)
    plt.imshow(self.fg_sigma)
    plt.colorbar()
    plt.title("foreground standard deviation")

    plt.subplot(232)
    plt.imshow(self.bg_mu,cmap=cm.gray,vmin=0,vmax=255)
    plt.colorbar()
    plt.title("background mean")

    plt.subplot(234)
    plt.imshow(self.bg_sigma)
    plt.colorbar()
    plt.title("background standard deviation")
    
    plt.subplot(235)
    plt.imshow(self.always_bg_mask_frac)
    plt.colorbar()
    plt.title("always bg mask frac")

    plt.figure()
    sz = (self.nr,self.nc,1)
    for j in range(4):
        # choose a frame
        self.frame = self.framessample[j]
        # read frame
        self.im,self.timestamp = self.movie.get_frame(self.frame)
        # read target positions
        self.trxcurr = self.trx.get_frame(self.frame)
        # create isfore mask
        self.isfore_mask()
        # detect objects
        self.obj_detect()
        [y,x] = num.where(self.isobj)
        plt.subplot(1,4,j+1)
        # red = image
        # green = 0
        # blue = isfore
        tmp = num.concatenate((self.im.reshape(sz),num.zeros(sz,dtype=num.uint8),255*self.isfore.astype(num.uint8).reshape(sz)),2)
        plt.imshow(tmp)
        plt.hold('on')
        plt.plot(x,y,'g.')
        plt.title("frame %d"%self.frame)
        plt.axis('image')


    plt.figure()
    plt.plot(self.obj_detection_dist_centers,self.obj_detection_dist_frac_bg,'k.-',
             self.obj_detection_dist_centers,self.obj_detection_dist_frac_fg,'r.-')
    plt.legend(('bg','fg'))
    plt.xlabel('Dist to obj detection')
    plt.ylabel('P(dist | label)')

    plt.show()

    return self

def test2():
    """
    Apply LoG filtering to a sample frame and show the results
    """

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
    
if __name__ == "__main__":
    test1()
