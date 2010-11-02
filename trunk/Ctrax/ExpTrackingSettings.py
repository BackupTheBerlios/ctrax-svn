import numpy as num
import scipy.stats as stats
from params import params
import annfiles as annot
import matplotlib.pyplot as plt
import os.path
import pickle
import getopt
import sys

if params.interactive:
    import wx

DEBUG = True
DELTAPLOT = .05

class ExpTrackingSettings:

    def __init__( self, annnames=None,
                  picklefile=None, 
                  nframessample = 100,
                  min_major_percentile = 1.,
                  max_major_percentile = 99.,
                  min_minor_percentile = 1.,
                  max_minor_percentile = 99.,
                  min_area_percentile = 1.,
                  max_area_percentile = 99.,
                  min_ecc_percentile = 1.,
                  max_ecc_percentile = 99.,
                  jump_distance_delta = .1
                  ):
        
        if picklefile is not None:
            if params.interactive:
                wx.Yield()
                wx.BeginBusyCursor()
            else:
                print 'Loading pickled tracking settings models'
            fid = open(picklefile,'r')
            model = pickle.load(fid)
            fid.close()

            for key,val in model.iteritems():
                setattr(self,key,val)

            if params.interactive:
                wx.EndBusyCursor()
        else:
        
            self.annnames = annnames
            self.nmovies = len(annnames)
    
            self.nframessample = nframessample
            self.min_major_percentile = min_major_percentile
            self.max_major_percentile = max_major_percentile
            self.min_minor_percentile = min_minor_percentile
            self.max_minor_percentile = max_minor_percentile
            self.min_area_percentile = min_area_percentile
            self.max_area_percentile = max_area_percentile
            self.min_ecc_percentile = min_ecc_percentile
            self.max_ecc_percentile = max_ecc_percentile
            self.jump_distance_delta = jump_distance_delta
            
            # initialize models
            self.init_shape_models()
            self.init_motion_models()

    def init_shape_models(self):

        self.majors = []
        self.minors = []
        self.areas = []
        self.eccs = []
        self.movieis = []
        self.flys = []
        self.frames = []

    def init_motion_models(self):

        self.xcurrs = []
        self.ycurrs = []
        self.anglecurrs = []
        self.xprevs = []
        self.yprevs = []
        self.angleprevs = []
        self.dxs = []
        self.dys = []
        self.dangles = []
        self.is_motion_data = []

    def update_shape_models(self):

        for ell in self.trxcurr.itervalues():
            self.majors.append(ell.size.height*4.)
            self.minors.append(ell.size.width*4.)
            if not hasattr(ell,'area') or ell.area < 0:
                ell.compute_area()
            self.areas.append(ell.area)
            self.eccs.append(ell.size.width / ell.size.height)
            self.movieis.append(self.moviei)
            self.flys.append(ell.identity)
            self.frames.append(self.frame)

    def update_motion_models(self):

        for id,ellcurr in self.trxcurr.iteritems():
            if self.trxprev.hasItem(id) and self.trxprevprev.hasItem(id):
                ellprev = self.trxprev[id]
                ellprevprev = self.trxprevprev[id]
                dx = ellprev.center.x - ellprevprev.center.x
                dy = ellprev.center.y - ellprevprev.center.y
                dangle = ((ellprev.angle - ellprevprev.angle + num.pi/2.) \
                              % (num.pi)) - (num.pi/2.)
                anglepred = ellprev.angle + dangle
                anglecurr = (ellcurr.angle - anglepred + num.pi/2.) % num.pi + anglepred - num.pi/2
                self.xcurrs.append(ellcurr.center.x)
                self.ycurrs.append(ellcurr.center.y)
                self.xprevs.append(ellprev.center.x)
                self.yprevs.append(ellprev.center.y)
                self.anglecurrs.append(anglecurr)
                self.angleprevs.append(ellprev.angle)
                self.dxs.append(dx)
                self.dys.append(dy)
                self.dangles.append(dangle)
                self.is_motion_data.append(True)

            else:

                self.xcurrs.append(num.nan)
                self.ycurrs.append(num.nan)
                self.xprevs.append(num.nan)
                self.yprevs.append(num.nan)
                self.anglecurrs.append(num.nan)
                self.angleprevs.append(num.nan)
                self.dxs.append(num.nan)
                self.dys.append(num.nan)
                self.dangles.append(num.nan)
                self.is_motion_data.append(False)
    
    def compute_models_permovie( self, i ):
        """
        compute_models_permovie( i )
        
        For annnames[i], this function samples frames from throughout
        the movie and computes the data structures necessary for 
        estimating the tracking settings.
        """

        self.moviei = i
        # open movie
        #self.movie = movies.Movie( self.movienames[i], params.interactive )

        # background model
        #self.bg_imgs = bg.BackgroundCalculator(self.movie)

        # open annotation
        self.trx = annot.AnnotationFile( self.annnames[i], doreadbgmodel=False, 
                                         doreadtrx=True )


        # choose frames to learn from: for now, assume all frames are tracked
        self.firstframe = self.trx.firstframetracked+2
        self.lastframe = self.trx.lastframetracked

        self.framessample = num.unique(num.round(num.linspace(self.firstframe+1,
                                                              self.lastframe,
                                                              self.nframessample)).astype(int))

        if DEBUG: print "Collecting data for movie %d, annfile = %s..."%(self.moviei,self.annnames[self.moviei])

        for j in range(len(self.framessample)):

            # read in the sampled frame
            self.frame = self.framessample[j]
            self.trxcurr = self.trx.get_frame(self.frame)
            self.trxprev = self.trx.get_frame(self.frame-1)
            self.trxprevprev = self.trx.get_frame(self.frame-2)

            if DEBUG: print "Collecting data for frame %d"%self.frame

            # update all the observation models
            self.update_shape_models()
            self.update_motion_models()

        #self.movie.close()
        #self.movie = None
        self.trx.close()
        self.trx = None
        #self.bg_imgs = None
        
    def est_settings(self):
        """
        est_settings()

        """

        for i in range(self.nmovies):

            if DEBUG: print "Computing per-movie model for movie %d"%i
            self.compute_models_permovie(i)

        self.est_shape_parameters()
        self.est_motion_parameters()

    def est_shape_parameters(self):

        self.majors = num.array(self.majors)
        self.minors = num.array(self.minors)
        self.areas = num.array(self.areas)
        self.eccs = num.array(self.eccs)

        self.min_major = genprctile(self.majors,self.min_major_percentile)
        self.max_major = genprctile(self.majors,self.max_major_percentile)
        self.min_minor = genprctile(self.minors,self.min_minor_percentile)
        self.max_minor = genprctile(self.minors,self.max_minor_percentile)
        self.min_area = genprctile(self.areas,self.min_area_percentile)
        self.max_area = genprctile(self.areas,self.max_area_percentile)
        self.min_ecc = genprctile(self.eccs,self.min_ecc_percentile)
        self.max_ecc = genprctile(self.eccs,self.max_ecc_percentile)
        self.mean_major = num.mean(self.majors)
        self.mean_minor = num.mean(self.minors)
        self.mean_area = num.mean(self.areas)
        self.mean_ecc = num.mean(self.eccs)

    def est_motion_parameters(self):

        self.xcurrs = num.array(self.xcurrs)
        self.ycurrs = num.array(self.ycurrs)
        self.anglecurrs = num.array(self.anglecurrs)
        self.xprevs = num.array(self.xprevs)
        self.yprevs = num.array(self.yprevs)
        self.angleprevs = num.array(self.angleprevs)
        self.dxs = num.array(self.dxs)
        self.dys = num.array(self.dys)
        self.dangles = num.array(self.dangles)
        self.is_motion_data = num.array(self.is_motion_data,dtype=bool)

        # min( [ (xcurrs-xprevs) - alpha*dxs ].T * (xcurrs-xprevs) - alpha*dxs ] + \
        #      [ (ycurrs-yprevs) - alpha*dys ].T * (ycurrs-yprevs) - alpha*dys ] )
        # =>
        # ((xcurrs-xprevs) - alpha*dxs).T*dxs + 
        # ((ycurrs-yprevs) - alpha*dxs).T*dxs = 0
        # (xcurrs-xprevs).T * dxs + (ycurrs-yprevs).T * dys = alpha*(dxs.T * dxs + dys.T * dy)
        # alpha = [ (xcurrs-xprevs).T * dxs + (ycurrs-yprevs).T * dys ] / (dxs.T * dxs + dys.T * dy)
        alpha = ( num.sum( (self.xcurrs[self.is_motion_data]-self.xprevs[self.is_motion_data])*self.dxs[self.is_motion_data] ) + \
                  num.sum( (self.ycurrs[self.is_motion_data]-self.yprevs[self.is_motion_data])*self.dys[self.is_motion_data] ) ) / \
                ( num.sum( self.dxs[self.is_motion_data]**2. ) + num.sum( self.dys[self.is_motion_data]**2. ) )
        alpha = max(alpha,0.)
        alpha = min(alpha,1.)
        self.center_dampen = 1. - alpha

        alpha = num.sum( self.anglecurrs[self.is_motion_data]*self.dangles[self.is_motion_data] ) / \
                num.sum( self.dangles[self.is_motion_data]**2. )
        alpha = max(alpha,0.)
        alpha = min(alpha,1.)
        self.angle_dampen = 1. - alpha

        # choose the weight of angle error
        self.xpreds = self.xprevs + self.dxs*(1.-self.center_dampen)
        self.ypreds = self.yprevs + self.dys*(1.-self.center_dampen)
        self.anglepreds = self.angleprevs + self.dangles*(1.-self.angle_dampen)
        self.center_err2s = (self.xpreds-self.xcurrs)**2 + (self.ypreds-self.ycurrs)**2.
        self.center_err2 = num.mean(self.center_err2s[self.is_motion_data])
        self.angle_err2s = ((self.anglepreds-self.anglecurrs+num.pi/2.)%num.pi - num.pi/2.)**2.
        self.angle_err2 = num.mean(self.angle_err2s[self.is_motion_data])
        self.angle_weight = self.center_err2 / self.angle_err2

        # choose the maximum jump distance
        self.dists = num.sqrt(self.center_err2s +self.angle_weight*self.angle_err2s)
        self.max_jump_distance = num.max(self.dists[self.is_motion_data])*(1.+self.jump_distance_delta)

    def save(self,outputFileName):
        
        fid = open(outputFileName,"w")
        
        out = dict()
        if hasattr(self,'nframessample'):
            out['nframessample'] = self.nframessample
        if hasattr(self,'min_major_percentile'):
            out['min_major_percentile'] = self.min_major_percentile
        if hasattr(self,'max_major_percentile'):
            out['max_major_percentile'] = self.max_major_percentile
        if hasattr(self,'min_minor_percentile'):
            out['min_minor_percentile'] = self.min_minor_percentile
        if hasattr(self,'max_minor_percentile'):
            out['max_minor_percentile'] = self.max_minor_percentile
        if hasattr(self,'min_area_percentile'):
            out['min_area_percentile'] = self.min_area_percentile
        if hasattr(self,'max_area_percentile'):
            out['max_area_percentile'] = self.max_area_percentile
        if hasattr(self,'min_ecc_percentile'):
            out['min_ecc_percentile'] = self.min_ecc_percentile
        if hasattr(self,'max_ecc_percentile'):
            out['max_ecc_percentile'] = self.max_ecc_percentile
        if hasattr(self,'jump_distance_delta'):
            out['jump_distance_delta'] = self.jump_distance_delta
        if hasattr(self,'annnames'):
            out['annnames'] = self.annnames

        if hasattr(self,'majors'):
            out['majors'] = self.majors
        if hasattr(self,'minors'):
            out['minors'] = self.minors
        if hasattr(self,'areas'):
            out['areas'] = self.areas
        if hasattr(self,'eccs'):
            out['eccs'] = self.eccs
        if hasattr(self,'movieis'):
            out['movieis'] = self.movieis
        if hasattr(self,'flys'):
            out['flys'] = self.flys
        if hasattr(self,'frames'):
            out['frames'] = self.frames
        if hasattr(self,'xcurrs'):
            out['xcurrs'] = self.xcurrs
        if hasattr(self,'ycurrs'):
            out['ycurrs'] = self.ycurrs
        if hasattr(self,'anglecurrs'):
            out['anglecurrs'] = self.anglecurrs
        if hasattr(self,'xprevs'):
            out['xprevs'] = self.xprevs
        if hasattr(self,'yprevs'):
            out['yprevs'] = self.yprevs
        if hasattr(self,'angleprevs'):
            out['angleprevs'] = self.angleprevs
        if hasattr(self,'dxs'):
            out['dxs'] = self.dxs
        if hasattr(self,'dys'):
            out['dys'] = self.dys
        if hasattr(self,'dangles'):
            out['dangles'] = self.dangles
        if hasattr(self,'is_motion_data'):
            out['is_motion_data'] = self.is_motion_data
        if hasattr(self,'min_major'):
            out['min_major'] = self.min_major
        if hasattr(self,'max_major'):
            out['max_major'] = self.max_major
        if hasattr(self,'min_minor'):
            out['min_minor'] = self.min_minor
        if hasattr(self,'max_minor'):
            out['max_minor'] = self.max_minor
        if hasattr(self,'min_area'):
            out['min_area'] = self.min_area
        if hasattr(self,'max_area'):
            out['max_area'] = self.max_area
        if hasattr(self,'min_ecc'):
            out['min_ecc'] = self.min_ecc
        if hasattr(self,'max_ecc'):
            out['max_ecc'] = self.max_ecc
        if hasattr(self,'mean_major'):
            out['mean_major'] = self.mean_major
        if hasattr(self,'mean_minor'):
            out['mean_minor'] = self.mean_minor
        if hasattr(self,'mean_area'):
            out['mean_area'] = self.mean_area
        if hasattr(self,'mean_ecc'):
            out['mean_ecc'] = self.mean_ecc
        if hasattr(self,'center_dampen'):
            out['center_dampen'] = self.center_dampen
        if hasattr(self,'angle_dampen'):
            out['angle_dampen'] = self.angle_dampen
        if hasattr(self,'xpreds'):
            out['xpreds'] = self.xpreds
        if hasattr(self,'ypreds'):
            out['ypreds'] = self.ypreds
        if hasattr(self,'anglepreds'):
            out['anglepreds'] = self.anglepreds
        if hasattr(self,'center_err2s'):
            out['center_err2s'] = self.center_err2s
        if hasattr(self,'angle_err2s'):
            out['angle_err2s'] = self.angle_err2s
        if hasattr(self,'center_err2'):
            out['center_err2'] = self.center_err2
        if hasattr(self,'angle_err2'):
            out['angle_err2'] = self.angle_err2
        if hasattr(self,'angle_weight'):
            out['angle_weight'] = self.angle_weight
        if hasattr(self,'dists'):
            out['dists'] = self.dists
        if hasattr(self,'max_jump_distance'):
            out['max_jump_distance'] = self.max_jump_distance

        pickle.dump(out,fid)
        fid.close()

    def show(self):
        
        self.show_shape()
        self.show_motion()
        plt.show()

    def show_motion(self):
        
        nbins = 100

        self.counts_centererr,self.edges_centererr = num.histogram(num.sqrt(self.center_err2s[self.is_motion_data]),nbins)
        self.counts_centererr = self.counts_centererr.astype(float)
        self.frac_centererr = self.counts_centererr / sum(self.counts_centererr)
        self.centers_centererr = (self.edges_centererr[:-1]+self.edges_centererr[1:])/2.

        self.counts_angleerr,self.edges_angleerr = num.histogram(num.sqrt(self.angle_err2s[self.is_motion_data]),nbins)
        self.counts_angleerr = self.counts_angleerr.astype(float)
        self.frac_angleerr = self.counts_angleerr / sum(self.counts_angleerr)
        self.centers_angleerr = (self.edges_angleerr[:-1]+self.edges_angleerr[1:])/2.

        self.counts_dist,self.edges_dist = num.histogram(self.dists[self.is_motion_data],nbins)
        self.counts_dist = self.counts_dist.astype(float)
        self.frac_dist = self.counts_dist / num.sum(self.counts_dist)
        self.centers_dist = (self.edges_dist[:-1]+self.edges_dist[1:])/2.

        plt.figure()

        plt.subplot(131)        
        plt.plot(self.centers_centererr,self.frac_centererr,'k.-')
        maxy = num.max(self.frac_centererr)
        ax = num.array([self.edges_centererr[0],self.edges_centererr[-1],-maxy*DELTAPLOT,maxy*(1.+DELTAPLOT)])
        plt.axis(ax)
        plt.title('center dist')

        plt.subplot(132)        
        plt.plot(self.centers_angleerr,self.frac_angleerr,'k.-')
        maxy = num.max(self.frac_angleerr)
        ax = num.array([self.edges_angleerr[0],self.edges_angleerr[-1],-maxy*DELTAPLOT,maxy*(1.+DELTAPLOT)])
        plt.axis(ax)
        plt.title('angle dist')
        
        plt.subplot(133)        
        plt.plot(self.centers_dist,self.frac_dist,'k.-')
        plt.hold(True)
        maxy = num.max(self.frac_dist)
        ax = num.array([self.edges_dist[0],max(self.max_jump_distance*(1.+DELTAPLOT),self.edges_dist[-1]),-maxy*DELTAPLOT,maxy*(1.+DELTAPLOT)])
        plt.plot(self.max_jump_distance+num.zeros(2),ax[[2,3]],'r-')
        plt.axis(ax)
        plt.title('dist')

        print "angle weight = " + str(self.angle_weight)
        print "center dampen = " + str(self.center_dampen)
        print "angle dampen = " + str(self.angle_dampen)
        print "max jump dist = " + str(self.max_jump_distance)

    def show_shape(self):

        nbins = 100

        # show shape
        
        # histogram of area axis lengths
        self.counts_area,self.edges_area = num.histogram(self.areas,nbins)
        self.counts_area = self.counts_area.astype(float)
        self.frac_area = self.counts_area / num.sum(self.counts_area)
        self.centers_area = (self.edges_area[:-1]+self.edges_area[1:])/2.

        # plot histogram
        plt.figure()

        plt.subplot(221)
        plt.plot(self.centers_area,self.frac_area,'k.-')
        plt.hold(True)
        ax = get_axis(self.edges_area,self.frac_area,self.min_area,self.max_area)
        plt.plot(self.min_area+num.zeros(2),ax[[2,3]],'r-')
        plt.plot(self.mean_area+num.zeros(2),ax[[2,3]],'r-')
        plt.plot(self.max_area+num.zeros(2),ax[[2,3]],'r-')
        plt.axis(ax)
        plt.title('area')

        # histogram of major axis lengths
        self.counts_major,self.edges_major = num.histogram(self.majors,nbins)
        self.counts_major = self.counts_major.astype(float)
        self.frac_major = self.counts_major / num.sum(self.counts_major)
        self.centers_major = (self.edges_major[:-1]+self.edges_major[1:])/2.

        # plot histogram
        plt.subplot(222)
        plt.plot(self.centers_major,self.frac_major,'k.-')
        plt.hold(True)
        ax = get_axis(self.edges_major,self.frac_major,self.min_major,self.max_major)
        plt.plot(self.min_major+num.zeros(2),ax[[2,3]],'r-')
        plt.plot(self.mean_major+num.zeros(2),ax[[2,3]],'r-')
        plt.plot(self.max_major+num.zeros(2),ax[[2,3]],'r-')
        plt.axis(ax)
        plt.title('major')

        # histogram of minor axis lengths
        self.counts_minor,self.edges_minor = num.histogram(self.minors,nbins)
        self.counts_minor = self.counts_minor.astype(float)
        self.frac_minor = self.counts_minor / num.sum(self.counts_minor)
        self.centers_minor = (self.edges_minor[:-1]+self.edges_minor[1:])/2.

        # plot histogram
        plt.subplot(223)
        plt.plot(self.centers_minor,self.frac_minor,'k.-')
        plt.hold(True)
        ax = get_axis(self.edges_minor,self.frac_minor,self.min_minor,self.max_minor)
        plt.plot(self.min_minor+num.zeros(2),ax[[2,3]],'r-')
        plt.plot(self.mean_minor+num.zeros(2),ax[[2,3]],'r-')
        plt.plot(self.max_minor+num.zeros(2),ax[[2,3]],'r-')
        plt.axis(ax)
        plt.title('minor')

        # histogram of ecc axis lengths
        self.counts_ecc,self.edges_ecc = num.histogram(self.eccs,nbins)
        self.counts_ecc = self.counts_ecc.astype(float)
        self.frac_ecc = self.counts_ecc / num.sum(self.counts_ecc)
        self.centers_ecc = (self.edges_ecc[:-1]+self.edges_ecc[1:])/2.

        # plot histogram
        plt.subplot(224)
        plt.plot(self.centers_ecc,self.frac_ecc,'k.-')
        plt.hold(True)
        ax = get_axis(self.edges_ecc,self.frac_ecc,self.min_ecc,self.max_ecc)
        plt.plot(self.min_ecc+num.zeros(2),ax[[2,3]],'r-')
        plt.plot(self.mean_ecc+num.zeros(2),ax[[2,3]],'r-')
        plt.plot(self.max_ecc+num.zeros(2),ax[[2,3]],'r-')
        plt.axis(ax)
        plt.title('ecc')

        print "min area = " + str(self.min_area)
        print "mean area = " + str(self.mean_area)
        print "max area = " + str(self.max_area)
        print "min major = " + str(self.min_major)
        print "mean major = " + str(self.mean_major)
        print "max major = " + str(self.max_major)
        print "min minor = " + str(self.min_minor)
        print "mean minor = " + str(self.mean_minor)
        print "max minor = " + str(self.max_minor)
        print "min ecc = " + str(self.min_ecc)
        print "mean ecc = " + str(self.mean_ecc)
        print "max ecc = " + str(self.max_ecc)

def main():
    """
    main()

    Compute experiment shape and motion models.
    """

    shortopts = "f:p:a:o:"
    longopts = ["filename=","params=","ann=","output="]
    try:
        opts, args = getopt.getopt(sys.argv[1:], shortopts, longopts)
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        print "ExpTrackingSettings options:"
        print "    -f,--filename <name of file containing list of experiment directories>"
        print "      Default: expdirs.txt"
        print "    -p,--params <name of file containing parameters>"
        print "      Default: <empty string>"
        print "    -a,--ann <name of annotation file within experiment directories>"
        print "      Default: movie.ufmf.ann"
        print "    -o,--output <name of file to output results to>"
        print "      Default: ExpTrackingSettingsResults.pickle"
        sys.exit(2)
        
    expdirsFileName = 'expdirs.txt'
    paramsFileName = ''
    annFileStr = 'movie.ufmf.ann'
    outputFileName = 'ExpTrackingSettingsResults.pickle'
    
    for o,a in opts:
        if o in ("-f","--filename"):
            expdirsFileName = a
            print "expdirsFileName = " + a
        elif o in ("-p","--params"):
            paramsFileName = a
            print "paramsFileName = " + a
        elif o in ("-a","--ann"):
            annFileStr = a
            print "annFileStr = " + a
        elif o in ("-o","--output"):
            outputFileName = a
            print "outputFileName = " + a
        else:
            assert False, "unhandled option"

    # read in the experiment directories
    fid = open(expdirsFileName,"r")
    expdirs = []
    annnames = []
    for l in fid:
        expdir = l.strip()
        if not os.path.exists(expdir):
            print "Experiment directory %s does not exist. Skipping."%expdir
            continue
        annname = os.path.join(expdir,annFileStr)
        if not os.path.exists(annname):
            print "Annotation file %s does not exist. Skipping experiment %s."%(annname,expdir)
            continue
        expdirs.append(expdir)
        annnames.append(annname)
    fid.close()
        
    if paramsFileName != '':
        params = read_params(paramsFileName)
        
    model = ExpTrackingSettings(annnames,
                                nframessample = params['nframessample'],
                                min_major_percentile = params['min_major_percentile'],
                                max_major_percentile = params['max_major_percentile'],
                                min_minor_percentile = params['min_minor_percentile'],
                                max_minor_percentile = params['max_minor_percentile'],
                                min_area_percentile = params['min_area_percentile'],
                                max_area_percentile = params['max_area_percentile'],
                                min_ecc_percentile = params['min_ecc_percentile'],
                                max_ecc_percentile = params['max_ecc_percentile'],
                                jump_distance_delta = params['jump_distance_delta'])
    model.est_settings()
    model.save(outputFileName)
    model.show()

    return model
    
def genprctile(x,prct):
    if prct < 0 or prct > 100:
        maxx = num.max(x)
        minx = num.min(x)
        dx = maxx - minx
        if prct < 0:
            p = minx - (-prct/100.)*dx
        else:
            p = maxx + ((prct-100.)/100)*dx
    else:
        p = stats.scoreatpercentile(x,prct)
    return p

def read_params(paramsFileName):

    # default parameters
    params = dict()
    params['nframessample'] = 100
    params['min_major_percentile'] = 1.
    params['max_major_percentile'] = 99.
    params['min_minor_percentile'] = 1.
    params['max_minor_percentile'] = 99.
    params['min_area_percentile'] = 1.
    params['max_area_percentile'] = 99.
    params['min_ecc_percentile'] = 1.
    params['max_ecc_percentile'] = 99.
    params['jump_distance_delta'] = .1

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
        if l[0] in params.keys():
            params[l[0]] = float(l[1])
    fid.close()

    params['nframessample'] = int(params['nframessample'])

    return params

def get_axis(edges,frac,minv,maxv):
    ax = num.zeros(4)
    dv = edges[-1] - edges[0]
    maxfrac = num.max(frac)
    ax[0] = min(edges[0],minv-dv*DELTAPLOT)
    ax[1] = max(edges[-1],maxv+dv*DELTAPLOT)
    ax[2] = 0-maxfrac*DELTAPLOT
    ax[3] = maxfrac*(1.+DELTAPLOT)
    return ax


if __name__ == "__main__":
    params.interactive = False
    main()
