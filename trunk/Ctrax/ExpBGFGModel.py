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
import getopt
import sys
import os.path
import pickle
import scipy.io

if params.interactive:
    import wx

DEBUG = True

if DEBUG:
    import pdb
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm


class ExpBGFGModel:

    def __init__( self, movienames=None, annnames=None, picklefile=None,
                  LoG_hsize=21, 
                  LoG_sigma=5.,
                  LoG_nonmaxsup_sz=5, 
                  LoG_thresh=.5,
                  difftype='darkfglightbg',
                  obj_detection_dist_nbins=100,
                  isback_dilation_hsize=10,
                  fg_min_sigma = 1.,
                  bg_min_sigma = 1.,
                  min_always_bg_mask_frac = 1.,
                  min_frac_cc_covered = .01,
                  min_area_cc_covered = 1,
                  backsub_n_bg_std_thresh = 10.,
                  backsub_n_bg_std_thresh_low = 1.,
                  backsub_do_use_morphology = False,
                  backsub_opening_radius = 0,
                  backsub_closing_radius = 0
                  ):
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
        to compute isback [10]. 
        fg_min_sigma: Minimum standard deviation for the foreground
        pixel appearance [1]. 
        bg_min_sigma: Minimum standard deviation for the background 
        pixel appearance [1]. 
        min_frac_cc_covered: Minimum fraction of a foreground connected component that
        must be covered by the tracked ellipse in order for the component to be added
        as foreground [.01].
        min_area_cc_covered: Minimum number of pixels of a connected component that
        must be covered by the tracked ellipse in order for the component to be added
        as foreground [1].
        backsub_n_bg_std_thresh: High threshold for background subtraction for the 
        sampled frame [10].
        backsub_n_bg_std_thresh_low: Low threshold for background subtraction for the
        sampled frame [1].
        backsub_do_use_morphology: Whether to do morphology after background subtraction
        [False].
        backsub_opening_radius: Radius of opening structural element after background
        subtraction [0].
        backsub_closing_radius: Radius of closing structural element after background
        subtraction [0].         
        """
        
        if picklefile is not None:
            if params.interactive:
                wx.Yield()
                wx.BeginBusyCursor()
            else:
                print 'Loading pickled global bg/fg models'
            fid = open(picklefile,'r')
            model = pickle.load(fid)
            fid.close()
            model['obj_detection_dist_centers'] = \
                (model['obj_detection_dist_edges'][1:] + \
                 model['obj_detection_dist_edges'][:-1])/2.

            for key,val in model.iteritems():
                setattr(self,key,val)
            (self.nr,self.nc) = self.bg_mu.shape

            if params.interactive:
                wx.EndBusyCursor()
        else:

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
            self.min_frac_cc_covered = min_frac_cc_covered
            self.min_area_cc_covered = min_area_cc_covered
            self.backsub_n_bg_std_thresh = backsub_n_bg_std_thresh
            self.backsub_n_bg_std_thresh_low = backsub_n_bg_std_thresh_low
            self.backsub_do_use_morphology = backsub_do_use_morphology
            self.backsub_opening_radius = backsub_opening_radius
            self.backsub_closing_radius = backsub_closing_radius

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

        # log( sqrt(2*pi) * sigma ) = .5*log(2*pi) + log(sigma)
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


        self.setup_backsub()

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

        self.movie.close()
        self.movie = None
        self.trx.close()
        self.trx = None
        self.bg_imgs = None
                                
    def compute_log_lik_appearance_given_fore(self,r0=0,r1=num.inf):
        """
        compute_log_lik_appearance_given_fore()
        
        computes the log likelihood of each pixel appearance in
        self.im given that the pixel is foreground. 
        """
        r1 = min(r1,self.nr)
        d2 = (self.im[r0:r1] - self.fg_mu[r0:r1])**2 / (2*self.fg_sigma[r0:r1]**2)
        self.log_lik_appearance_given_fore = -d2 - self.fg_log_Z[r0:r1]

    def compute_log_lik_appearance_given_back(self,r0=0,r1=num.inf):
        """
        compute_log_lik_appearance_given_back()
        
        computes the log likelihood of each pixel appearance in
        self.im given that the pixel is background. 
        """
        
        r1 = min(r1,self.nr)

        d2 = (self.im[r0:r1] - self.bg_mu[r0:r1])**2 / (2*self.bg_sigma[r0:r1]**2)
        self.log_lik_appearance_given_back = -d2 - self.bg_log_Z[r0:r1]

    def compute_log_lik_dist_obj(self,r0=0,r1=num.inf):
        """
        compute_log_lik_dist_obj()
        
        computes the log likelihood of the distance to the nearest object detection
        for each pixel in the image self.im given that the pixel is in foreground and
        given that the pixel is in background. 
        """
        
        r1 = min(r1,self.nr)
        
        # detect objects in the current image
        self.obj_detect(r0=r0,r1=r1)

        # compute distances to detections
        morph.distance_transform_edt(~self.isobj,
                                     distances=self.obj_detection_dist[r0:r1])
        nr = r1 - r0

        # find which bin each distance falls in
        idx = num.digitize(num.reshape(self.obj_detection_dist[r0:r1],(nr*self.nc,)),self.obj_detection_dist_edges)
        idx.shape = (nr,self.nc)

        # lookup probabilities for each bin idx
        self.log_lik_dist_obj_given_fore = num.log(self.obj_detection_dist_frac_fg[idx])
        self.log_lik_dist_obj_given_back = num.log(self.obj_detection_dist_frac_bg[idx])

    def compute_log_lik(self,r0=0,r1=num.inf):

        self.compute_log_lik_appearance_given_fore(r0=r0,r1=r1)
        self.compute_log_lik_appearance_given_back(r0=r0,r1=r1)
        self.compute_log_lik_dist_obj(r0=r0,r1=r1)
        self.log_lik_given_fore = self.log_lik_appearance_given_fore + \
            self.log_lik_dist_obj_given_fore
        self.log_lik_given_back = self.log_lik_appearance_given_back + \
            self.log_lik_dist_obj_given_back
            
    def compute_log_lik_ratio(self,r0=0,r1=num.inf):
        
        self.compute_log_lik(r0=r0,r1=r1)
        self.log_lik_ratio = self.log_lik_given_fore - self.log_lik_given_back
        return self.log_lik_ratio
            
    def setup_backsub(self):

        params.n_bg_std_thresh = self.backsub_n_bg_std_thresh
        params.n_bg_std_thresh_low = self.backsub_n_bg_std_thresh_low
        params.do_use_morphology = self.backsub_do_use_morphology
        params.opening_radius = self.backsub_opening_radius
        params.closing_radius = self.backsub_closing_radius
        if params.do_use_morphology:
            self.bg_imgs.opening_struct = self.bg_imgs.create_morph_struct(params.opening_radius)
            self.bg_imgs.closing_struct = self.bg_imgs.create_morph_struct(params.closing_radius)

    def backsub(self):

        return self.bg_imgs.sub_bg(self.frame)

    def isfore_mask(self):

        self.isfore = num.zeros((self.nr,self.nc),dtype=bool)
        idx = num.zeros((self.nr,self.nc),dtype=bool)
        #print 'backsub_n_bg_std_thresh = ' + str(params.n_bg_std_thresh)
        #print 'backsub_n_bg_std_thresh should be ' + str(self.backsub_n_bg_std_thresh)
        #print 'backsub_n_bg_std_thresh_low = ' + str(params.n_bg_std_thresh_low)
        #print 'backsub_do_use_morphology = ' + str(params.do_use_morphology)
        #print 'backsub_opening_radius = ' + str(params.opening_radius)
        #print 'backsub_closing_radius = ' + str(params.closing_radius)
        #print 'min_frac_cc_covered = ' + str(self.min_frac_cc_covered)
        #print 'min_area_cc_covered = ' + str(self.min_area_cc_covered)
        (self.dfore,self.bw,self.cc,self.ncc) = self.backsub()

        # area of each connected component
        areas = num.zeros(self.ncc)
        for i in range(self.ncc):
            areas[i] = num.alen(num.flatnonzero(self.cc==i+1))

        for ell in self.trxcurr.itervalues():

            # find pixels inside this ellipse
            (r1,r2,c1,c2) = est.getboundingboxtight(ell,(self.nr,self.nc))
            isfore1 = est.ellipsepixels(ell,num.array([r1,r2,c1,c2]))

            # which cc(s) does this correspond to?
            r,c = isfore1.nonzero()
            counts,edges = num.histogram(self.cc[int(r1)+r,int(c1)+c],num.arange(1,self.ncc+2))
            fraccovered = counts / areas
            iscovered = fraccovered >= self.min_frac_cc_covered
            mode = num.argmax(counts)
            if counts[mode] > self.min_area_cc_covered:
                iscovered[mode] = True
            #print "counts = " + str(counts)
            #print "areas = " + str(areas)
            #print "fraccovered = " + str(fraccovered)
            #print "mode = " + str(mode)

            # set isfore to True for these ccs
            for i in num.flatnonzero(iscovered):
                self.isfore[self.cc==i+1] = True

            #h0 = plt.subplot(131)
            #plt.imshow(self.im)
            #plt.hold(True)
            #plt.plot(ell.x,ell.y,'ro')
            #plt.axis('image')
            #plt.title('image')
            #plt.subplot(132,sharex=h0,sharey=h0)
            #plt.title('cc')
            #plt.imshow(self.cc)
            #plt.subplot(133,sharex=h0,sharey=h0)
            #isforecurr = num.zeros((self.nr,self.nc),dtype=bool)
            #for i in num.flatnonzero(iscovered):
            #    isforecurr[self.cc==i+1] = True
            #plt.imshow(isforecurr)
            #if num.any(isforecurr):
            #    rcurr,ccurr = num.where(isforecurr)
            #    r0curr = num.min(rcurr)-1
            #    r1curr = num.max(rcurr)+1
            #    c0curr = num.min(ccurr)-1
            #    c1curr = num.max(ccurr)+1
            #    plt.hold(True)
            #    plt.plot(num.array([c0curr,c0curr,c1curr,c1curr,c0curr]),
            #             num.array([r0curr,r1curr,r1curr,r0curr,r0curr]),'w-')
            #plt.title('isfore curr')
            #plt.show()


        # compute isback by dilating and flipping
        self.isback = morph.binary_dilation(self.isfore,self.dilation_struct) == False

        h0 = plt.subplot(151)
        plt.imshow(self.im)
        plt.axis('tight')
        plt.axis('image')
        plt.title('image')
        plt.subplot(152,sharex=h0,sharey=h0)
        plt.title('cc')
        plt.imshow(self.cc)
        plt.subplot(153,sharex=h0,sharey=h0)
        plt.imshow(self.isfore)
        plt.title('isfore')
        plt.subplot(154,sharex=h0,sharey=h0)
        plt.imshow(self.isback)
        plt.subplot(155,sharex=h0,sharey=h0)
        plt.imshow(num.logical_and(self.bw,self.isfore==False))
        plt.title('ignored fg pixels')
        plt.show()
            
    def isfore_mask_old(self):
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
                      
                      
    def obj_detect(self,r0=0,r1=num.inf):
        """
        obj_detect()

        Detect objects in the current frame. For now, this is just
        LoG filtering. 
        """

        r1 = min(r1,self.nr)

        self.LoGfilter(r0=r0,r1=r1)
        
        if hasattr(self,'always_bg_mask') and self.always_bg_mask is not None:
            self.isobj[self.always_bg_mask[r0:r1]] = False

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

        # histograms of distances to object detections for foreground and background:
        # initialize to 1 so that there is never a 0 probability
        self.obj_detection_dist_counts_fg = num.ones(self.obj_detection_dist_nbins)
        self.obj_detection_dist_counts_bg = num.ones(self.obj_detection_dist_nbins)

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
        

    def LoGfilter(self,r0=0,r1=num.inf):
        """
        LoGfilter()

        Apply LoG filtering to the current frame to detect blobs. 
        """
        
        r1 = min(r1,self.nr)
        
        # LoG filter
        if not hasattr(self,'LoG_fil_out'):
            # output of LoG filtering
            self.LoG_fil_out = num.zeros((self.nr,self.nc))
            
        if not hasattr(self,'LoG_fil_max'):
            # output of nonmaximal suppression
            self.LoG_fil_max = num.zeros((self.nr,self.nc))

        if not hasattr(self,'obj_detection_dist'):
            # output of dist transform
            self.obj_detection_dist = num.zeros((self.nr,self.nc))
            
        # pad borders
        off = max(self.LoG_hsize,self.LoG_nonmaxsup_sz) 
        r00 = max(0,r0-off)
        r11 = min(self.nr,r1+off)

        filters.correlate(self.im[r00:r11].astype(float),
                          self.LoG_fil,
                          output=self.LoG_fil_out[r00:r11],
                          mode='nearest')

        # depending on difftype, we will only look for positive, negative or both values
        if self.difftype == 'other':
            self.LoG_fil_out[r00:r11] = num.abs(self.LoG_fil_out[r00:r11])

        # non-maximal suppression + threshold
        filters.maximum_filter(self.LoG_fil_out[r00:r11],size=self.LoG_nonmaxsup_sz,
                               output=self.LoG_fil_max[r00:r11],mode='constant',cval=0)
        self.isobj = num.logical_and(self.LoG_fil_out[r0:r1] == self.LoG_fil_max[r0:r1],
                                     self.LoG_fil_out[r0:r1] >= self.LoG_thresh)
                                     
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
    
    def save_mat(self,matFileName):
        
        out = dict()
        if hasattr(self,'fg_mu'):
            out['fg_mu'] = self.fg_mu
        if hasattr(self,'fg_sigma'):
            out['fg_sigma'] = self.fg_sigma
        if hasattr(self,'fg_log_Z'):
            out['fg_log_Z'] = self.fg_log_Z
        if hasattr(self,'bg_mu'):
            out['bg_mu'] = self.bg_mu
        if hasattr(self,'bg_sigma'):
            out['bg_sigma'] = self.bg_sigma
        if hasattr(self,'bg_log_Z'):
            out['bg_log_Z'] = self.bg_log_Z
        if hasattr(self,'obj_detection_dist_frac_fg'):
            out['obj_detection_dist_frac_fg'] = self.obj_detection_dist_frac_fg
        if hasattr(self,'obj_detection_dist_frac_bg'):
            out['obj_detection_dist_frac_bg'] = self.obj_detection_dist_frac_bg
        if hasattr(self,'obj_detection_dist_edges'):
            out['obj_detection_dist_edges'] = self.obj_detection_dist_edges
        if hasattr(self,'obj_detection_dist_centers'):
            out['obj_detection_dist_centers'] = self.obj_detection_dist_centers
        if hasattr(self,'always_bg_mask_frac'):
            out['always_bg_mask_frac'] = self.always_bg_mask_frac
        if hasattr(self,'always_bg_mask'):
            out['always_bg_mask'] = self.always_bg_mask

        scipy.io.savemat( matfilename, out, oned_as='column' )

    def save(self,outputFileName):
        
        fid = open(outputFileName,"w")
        
        out = dict()
        if hasattr(self,'LoG_hsize'):
            out['LoG_hsize'] = self.LoG_hsize
        if hasattr(self,'LoG_sigma'):
            out['LoG_sigma'] = self.LoG_sigma
        if hasattr(self,'LoG_nonmaxsup_sz'):
            out['LoG_nonmaxsup_sz'] = self.LoG_nonmaxsup_sz
        if hasattr(self,'LoG_thresh'):
            out['LoG_thresh'] = self.LoG_thresh
        if hasattr(self,'difftype'):
            out['difftype'] = self.difftype
        if hasattr(self,'obj_detection_dist_nbins'):
            out['obj_detection_dist_nbins'] = self.obj_detection_dist_nbins
        if hasattr(self,'isback_dilation_hsize'):
            out['isback_dilation_hsize'] = self.isback_dilation_hsize
        if hasattr(self,'fg_min_sigma'):
            out['fg_min_sigma'] = self.fg_min_sigma
        if hasattr(self,'bg_min_sigma'):
            out['bg_min_sigma'] = self.bg_min_sigma
        if hasattr(self,'min_always_bg_mask_frac'):
            out['min_always_bg_mask_frac'] = self.min_always_bg_mask_frac
        if hasattr(self,'min_frac_cc_covered'):
            out['min_frac_cc_covered'] = self.min_frac_cc_covered
        if hasattr(self,'min_area_cc_covered'):
            out['min_area_cc_covered'] = self.min_area_cc_covered
        if hasattr(self,'backsub_n_bg_std_thresh'):
            out['backsub_n_bg_std_thresh'] = self.backsub_n_bg_std_thresh
        if hasattr(self,'backsub_n_bg_std_thresh_low'):
            out['backsub_n_bg_std_thresh_low'] = self.backsub_n_bg_std_thresh_low
        if hasattr(self,'backsub_do_use_morphology'):
            out['backsub_do_use_morphology'] = self.backsub_do_use_morphology
        if hasattr(self,'backsub_opening_radius'):
            out['backsub_opening_radius'] = self.backsub_opening_radius
        if hasattr(self,'backsub_closing_radius'):
            out['backsub_closing_radius'] = self.backsub_closing_radius
        if hasattr(params,'prior_nframessample'):
            out['prior_nframessample'] = params.prior_nframessample
        if hasattr(params,'prior_fg_nsamples_pool'):
            out['prior_fg_nsamples_pool'] = params.prior_fg_nsamples_pool
        if hasattr(params,'prior_bg_nsamples_pool'):
            out['prior_bg_nsamples_pool'] = params.prior_bg_nsamples_pool
        if hasattr(params,'prior_fg_pool_radius_factor'):
            out['prior_fg_pool_radius_factor'] = params.prior_fg_pool_radius_factor
        if hasattr(params,'prior_bg_pool_radius_factor'):
            out['prior_bg_pool_radius_factor'] = params.prior_bg_pool_radius_factor
        if hasattr(self,'movienames'):
            out['movienames'] = self.movienames
        if hasattr(self,'annnames'):
            out['annnames'] = self.annnames
        if hasattr(self,'fg_mu'):
            out['fg_mu'] = self.fg_mu
        if hasattr(self,'fg_sigma'):
            out['fg_sigma'] = self.fg_sigma
        if hasattr(self,'fg_log_Z'):
            out['fg_log_Z'] = self.fg_log_Z
        if hasattr(self,'bg_mu'):
            out['bg_mu'] = self.bg_mu
        if hasattr(self,'bg_sigma'):
            out['bg_sigma'] = self.bg_sigma
        if hasattr(self,'bg_log_Z'):
            out['bg_log_Z'] = self.bg_log_Z
        if hasattr(self,'obj_detection_dist_frac_fg'):
            out['obj_detection_dist_frac_fg'] = self.obj_detection_dist_frac_fg
        if hasattr(self,'obj_detection_dist_frac_bg'):
            out['obj_detection_dist_frac_bg'] = self.obj_detection_dist_frac_bg
        if hasattr(self,'obj_detection_dist_edges'):
            out['obj_detection_dist_edges'] = self.obj_detection_dist_edges
        if hasattr(self,'LoG_fil'):
            out['LoG_fil'] = self.LoG_fil
        if hasattr(self,'always_bg_mask_frac'):
            out['always_bg_mask_frac'] = self.always_bg_mask_frac
        if hasattr(self,'always_bg_mask'):
            out['always_bg_mask'] = self.always_bg_mask

        pickle.dump(out,fid)
        fid.close()


    def set_movie(self,moviename):

        self.movie = movies.Movie( moviename, params.interactive )
        self.nr = self.movie.get_height()
        self.nc = self.movie.get_width()


    def show(self):

        plt.subplot(231)
        plt.imshow(self.fg_mu,cmap=cm.gray,vmin=0,vmax=255)
        plt.colorbar()
        plt.title("foreground mean")
    
        plt.subplot(234)
        plt.imshow(self.fg_sigma)
        plt.colorbar()
        plt.title("foreground standard deviation")
    
        plt.subplot(232)
        plt.imshow(self.bg_mu,cmap=cm.gray,vmin=0,vmax=255)
        plt.colorbar()
        plt.title("background mean")
    
        plt.subplot(235)
        plt.imshow(self.bg_sigma)
        plt.colorbar()
        plt.title("background standard deviation")
        
        plt.subplot(233)
        plt.imshow(self.always_bg_mask_frac)
        plt.colorbar()
        plt.title("always bg mask frac")
    
        plt.subplot(236)
        plt.imshow(self.always_bg_mask,cmap=cm.gray)
        plt.colorbar()
        plt.title("always bg mask")
    
        plt.figure()
        plt.plot(self.obj_detection_dist_centers,self.obj_detection_dist_frac_bg,'k.-',
                 self.obj_detection_dist_centers,self.obj_detection_dist_frac_fg,'r.-')
        plt.legend(('bg','fg'))
        plt.xlabel('Dist to obj detection')
        plt.ylabel('P(dist | label)')
    
        plt.show()
        
    def showtest(self,moviename=None,moviei=0,nframessample=4,cmin=-num.infty,cmax=num.infty):
        
        if moviename is None:
            moviename = self.movienames[moviei]
            
        self.movie = movies.Movie( moviename, params.interactive )
        self.nr = self.movie.get_height()
        self.nc = self.movie.get_width()

        framessample = num.unique(num.round(num.linspace(0,self.movie.get_n_frames()-1,nframessample)).astype(int))
        
        nr = nframessample
        nc = 5
        
        plt.figure()
        vmin_noninf = num.zeros(nc)
        vmin_noninf[:] = num.infty
        vmax_noninf = num.zeros(nc)
        vmax_noninf[:] = -num.infty
        vmin = num.zeros(nc)
        vmin[:] = num.infty
        vmax = num.zeros(nc)
        vmax[:] = -num.infty
        
        hims = []
        for i in range(nframessample):
            
            self.frame = framessample[i]
            print 'Frame %d (%d/%d)'%(self.frame,i+1,nframessample)
            self.im,self.timestamp = self.movie.get_frame(self.frame)
            log_lik_ratio = self.compute_log_lik_ratio()
            
            j = 1
            h0 = plt.subplot(nr,nc,nc*i+j)
            tmp = self.log_lik_appearance_given_fore.copy()
            tmp[self.always_bg_mask] = num.nan
            him = []
            him.append(plt.imshow(tmp))
            plt.axis('off')
            plt.axis('image')
            plt.title('ll(intensity|fg), %d'%self.frame)
            (vmin,vmax,vmin_noninf,vmax_noninf) = update_clim(self.log_lik_appearance_given_fore,vmin,vmax,vmin_noninf,vmax_noninf,j)
            plt.colorbar()
            
            j += 1
            plt.subplot(nr,nc,nc*i+j,sharex=h0,sharey=h0)
            tmp = self.log_lik_appearance_given_back.copy()
            tmp[self.always_bg_mask] = num.nan
            him.append(plt.imshow(tmp))
            plt.axis('off')
            plt.title('ll(intensity|bg), %d'%self.frame)
            (vmin,vmax,vmin_noninf,vmax_noninf) = update_clim(self.log_lik_appearance_given_back,vmin,vmax,vmin_noninf,vmax_noninf,j)
            plt.colorbar()

            j += 1
            plt.subplot(nr,nc,nc*i+j,sharex=h0,sharey=h0)
            tmp = self.log_lik_dist_obj_given_fore.copy()
            tmp[self.always_bg_mask] = num.nan
            him.append(plt.imshow(tmp))
            plt.axis('off')
            plt.title('ll(distobj|fg), %d'%self.frame)
            (vmin,vmax,vmin_noninf,vmax_noninf) = update_clim(self.log_lik_dist_obj_given_fore,vmin,vmax,vmin_noninf,vmax_noninf,j)
            plt.colorbar()

            j += 1
            plt.subplot(nr,nc,nc*i+4,sharex=h0,sharey=h0)
            tmp = self.log_lik_dist_obj_given_back.copy()
            tmp[self.always_bg_mask] = num.nan
            him.append(plt.imshow(tmp))
            plt.axis('off')
            plt.title('ll(distobj|bg), %d'%self.frame)
            (vmin,vmax,vmin_noninf,vmax_noninf) = update_clim(self.log_lik_dist_obj_given_back,vmin,vmax,vmin_noninf,vmax_noninf,j)
            plt.colorbar()
            
            j += 1
            plt.subplot(nr,nc,nc*i+j,sharex=h0,sharey=h0)
            tmp = log_lik_ratio.copy()
            tmp[self.always_bg_mask] = num.nan
            him.append(plt.imshow(tmp))
            plt.axis('off')
            plt.title('ll(im|fg) - ll(im|bg), %d'%self.frame)
            (vmin,vmax,vmin_noninf,vmax_noninf) = update_clim(log_lik_ratio,vmin,vmax,vmin_noninf,vmax_noninf,j)
            plt.colorbar()
            
            hims.append(him)
            plt.draw()

        vmin = num.maximum(vmin,cmin)
        vmin_noninf = num.maximum(vmin_noninf,cmin)
        vmax = num.minimum(vmax,cmax)
        vmax_noninf = num.minimum(vmax_noninf,cmax)
        dv = vmax_noninf - vmin_noninf
        idx = vmin < vmin_noninf
        vmin_noninf[idx] = vmin_noninf[idx] - dv[idx]*.025
        idx = vmax > vmax_noninf
        vmax_noninf[idx] = vmax_noninf[idx] + dv[idx]*.025

        for j in range(nc):
            for i in range(nr):
                hims[i][j].set_clim(vmin_noninf[j],vmax_noninf[j])
                
        print "Done."
        
        plt.show()
            
def update_clim(x,vmin,vmax,vmin_noninf,vmax_noninf,j):
    
    vmin[j-1] = min(num.min(x),vmin[j-1])
    vmax[j-1] = max(num.max(x),vmax[j-1])
    vmin_noninf[j-1] = min(num.min(x[False == num.isinf(x)]),vmin_noninf[j-1])
    vmax_noninf[j-1] = max(num.max(x[False == num.isinf(x)]),vmax_noninf[j-1])

    return (vmin,vmax,vmin_noninf,vmax_noninf)
            
    
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

    plt.subplot(234)
    plt.imshow(self.fg_sigma)
    plt.colorbar()
    plt.title("foreground standard deviation")

    plt.subplot(232)
    plt.imshow(self.bg_mu,cmap=cm.gray,vmin=0,vmax=255)
    plt.colorbar()
    plt.title("background mean")

    plt.subplot(235)
    plt.imshow(self.bg_sigma)
    plt.colorbar()
    plt.title("background standard deviation")
    
    plt.subplot(233)
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

def main():
    """
    main()

    Compute experiment bg and fg models for a set of movies set
    at the start of the function, and plot the results. 
    """

    shortopts = "f:p:m:a:o:s:"
    longopts = ["filename=","params=","movie=","ann=","output=","mat="]
    try:
        opts, args = getopt.getopt(sys.argv[1:], shortopts, longopts)
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        print "ExpBGFGModel options:"
        print "    -f,--filename <name of file containing list of experiment directories>"
        print "      Default: expdirs.txt"
        print "    -p,--params <name of file containing parameters>"
        print "      Default: <empty string>"
        print "    -m,--movie <name of movie file within experiment directories>"
        print "      Default: movie.ufmf"
        print "    -a,--ann <name of annotation file within experiment directories>"
        print "      Default: movie.ufmf.ann"
        print "    -o,--output <name of file to output results to>"
        print "      Default: ExpBGFGModelResults.pickle"
        print "    --mat <name of mat file to output results to>"
        print "      Default: ExpBGFGModelResults.mat"
        sys.exit(2)
        
    expdirsFileName = 'expdirs.txt'
    paramsFileName = ''
    movieFileStr = 'movie.ufmf'
    annFileStr = 'movie.ufmf.ann'
    outputFileName = 'ExpBGFGModelResults.pickle'
    matFileName = 'ExpBGFGModelResults.mat'
    
    for o,a in opts:
        if o in ("-f","--filename"):
            expdirsFileName = a
            print "expdirsFileName = " + a
        elif o in ("-p","--params"):
            paramsFileName = a
            print "paramsFileName = " + a
        elif o in ("-m","--movie"):
            movieFileStr = a
            print "movieFileStr = " + a
        elif o in ("-a","--ann"):
            annFileStr = a
            print "annFileStr = " + a
        elif o in ("-o","--output"):
            outputFileName = a
            print "outputFileName = " + a
        elif o in ("-s","--mat"):
            matFileName = a
            print "matFileName = " + a
        else:
            assert False, "unhandled option"

    # read in the experiment directories
    fid = open(expdirsFileName,"r")
    expdirs = []
    movienames = []
    annnames = []
    for l in fid:
        expdir = l.strip()
        if not os.path.exists(expdir):
            print "Experiment directory %s does not exist. Skipping."%expdir
            continue
        moviename = os.path.join(expdir,movieFileStr)
        if not os.path.exists(moviename):
            print "Movie %s does not exist. Skipping experiment %s."%(moviename,expdir)
            continue
        annname = os.path.join(expdir,annFileStr)
        if not os.path.exists(annname):
            print "Annotation file %s does not exist. Skipping experiment %s."%(annname,expdir)
            continue
        expdirs.append(expdir)
        movienames.append(moviename)
        annnames.append(annname)
    fid.close()
        
    # default parameters
    LoG_hsize=21 
    LoG_sigma=5.
    LoG_nonmaxsup_sz=5 
    LoG_thresh=.5
    difftype='darkfglightbg'
    obj_detection_dist_nbins=100
    isback_dilation_hsize=10
    fg_min_sigma = 1.
    bg_min_sigma = 1.
    min_always_bg_mask_frac = 1. 
    min_frac_cc_covered = .01
    min_area_cc_covered = 1
    backsub_n_bg_std_thresh = 10.
    backsub_n_bg_std_thresh_low = 1.
    backsub_do_use_morphology = False
    backsub_opening_radius = 0
    backsub_closing_radius = 0

    if paramsFileName != '':
        # parse parameters file
        fid = open(paramsFileName,"r")
    
        for l in fid:
            l = l.strip()
            if l[0] == '#':
                continue
            l = l.split('=',1)
            if len(l) < 2:
                print "Skipping parameter line '%s'"%l
                continue
            l[0] = l[0].strip()
            l[1] = l[1].strip()
            if l[0] == 'LoG_hsize':
                LoG_hsize = float(l[1])
            elif l[0] == 'LoG_sigma':
                LoG_sigma = float(l[1])
            elif l[0] == 'LoG_nonmaxsup_sz':
                LoG_nonmaxsup_sz = float(l[1]) 
            elif l[0] == 'LoG_thresh':
                LoG_thresh = float(l[1])
            elif l[0] == 'difftype':
                difftype = l[1]
            elif l[0] == 'obj_detection_dist_nbins':
                obj_detection_dist_nbins = float(l[1])
            elif l[0] == 'isback_dilation_hsize':
                isback_dilation_hsize = float(l[1])
            elif l[0] == 'fg_min_sigma':
                fg_min_sigma = float(l[1])
            elif l[0] == 'bg_min_sigma':
                bg_min_sigma = float(l[1])
            elif l[0] == 'min_always_bg_mask_frac':
                min_always_bg_mask_frac = float(l[1]) 
            elif l[0] == 'min_frac_cc_covered':
                min_frac_cc_covered = float(l[1])
            elif l[0] == 'min_area_cc_covered':
                min_area_cc_covered = float(l[1])
            elif l[0] == 'backsub_n_bg_std_thresh':
                backsub_n_bg_std_thresh = float(l[1])
            elif l[0] == 'backsub_n_bg_std_thresh_low':
                backsub_n_bg_std_thresh_low = float(l[1])
            elif l[0] == 'backsub_do_use_morphology':
                backsub_do_use_morphology = float(l[1]) == 1
            elif l[0] == 'backsub_opening_radius':
                backsub_opening_radius = int(float(l[1]))
            elif l[0] == 'backsub_closing_radius':
                backsub_closing_radius = int(float(l[1]))
            elif l[0] == 'prior_nframessample':
                params.prior_nframessample = float(l[1])
            elif l[0] == 'prior_fg_nsamples_pool':
                params.prior_fg_nsamples_pool = float(l[1])
            elif l[0] == 'prior_bg_nsamples_pool':
                params.prior_bg_nsamples_pool = float(l[1])
            elif l[0] == 'prior_fg_pool_radius_factor':
                params.prior_fg_pool_radius_factor = float(l[1])
            elif l[0] == 'prior_bg_pool_radius_factor':
                params.prior_bg_pool_radius_factor = float(l[1])
        fid.close()

    model = ExpBGFGModel(movienames,annnames,
                        LoG_hsize=LoG_hsize, 
                        LoG_sigma=LoG_sigma,
                        LoG_nonmaxsup_sz=LoG_nonmaxsup_sz, 
                        LoG_thresh=LoG_thresh,
                        difftype=difftype,
                        obj_detection_dist_nbins=obj_detection_dist_nbins,
                        isback_dilation_hsize=isback_dilation_hsize,
                        fg_min_sigma=fg_min_sigma,
                        bg_min_sigma=bg_min_sigma,
                        min_always_bg_mask_frac=min_always_bg_mask_frac,
                        min_frac_cc_covered = min_frac_cc_covered,
                        min_area_cc_covered = min_area_cc_covered,
                        backsub_n_bg_std_thresh = backsub_n_bg_std_thresh,
                        backsub_n_bg_std_thresh_low = backsub_n_bg_std_thresh_low,
                        backsub_do_use_morphology = backsub_do_use_morphology,
                        backsub_opening_radius = backsub_opening_radius,
                        backsub_closing_radius = backsub_closing_radius
                        )
    
    model.est_marginals()
    model.save(outputFileName)
    model.save_mat(matFileName)
    model.show()

    return model
    
if __name__ == "__main__":
    params.interactive = False
    main()
