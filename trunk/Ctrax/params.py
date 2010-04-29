from version import DEBUG, __version__
import numpy as num
import copy

# about box
# forward compatibility for AboutDialogInfo, added in wx 2.7.1.1
# (currently 2.6.3.2 is packaged for Ubuntu Feisty)
try:
    #from wx import AboutDialogInfo
    raise ImportError # because the built-in versions have not been tested
except ImportError:
    from about import AboutDialogInfo

class ShapeParams:
    def __init__(self,major=0,minor=0,area=0,ecc=0):
        self.major = major
        self.minor = minor
        self.area = area
        self.ecc = ecc

    def copy(self):
        return ShapeParams(self.major,self.minor,self.area,self.ecc)

    def __print__(self):
        return 'major = %.2f, minor = %.2f, area = %.2f, ecc = %.2f'%(self.major,self.minor,self.area,self.ecc)

    def __repr__( self ): return self.__print__()

    def __str__( self ): return self.__print__()

    def __eq__(self,other):
        for i,j in self.__dict__.iteritems():
            if not hasattr(other,i):
                return False
            if not (j == other.__dict__[i]):
                return False
        for i,j in other.__dict__.iteritems():
            if not hasattr(self,i):
                return False
            if not (j == self.__dict__[i]):
                return False
        return True

def averageshape(shape1,shape2):
    shape3 = ShapeParams()
    shape3.major = (shape1.major+shape2.major)/2.
    shape3.minor = (shape1.minor+shape2.minor)/2.
    shape3.area = (shape1.area+shape2.area)/2.
    shape3.ecc = (shape1.ecc+shape2.ecc)/2.
    return shape3

class Grid:
    def __init__(self):
        self.X = None
        self.Y = None
        self.X2 = None
        self.Y2 = None
        self.XY = None
    def setsize(self,sz):
        [self.Y,self.X] = num.mgrid[0:sz[0],0:sz[1]]
        self.Y2 = self.Y**2
        self.X2 = self.X**2
        self.XY = self.X*self.Y
    def __eq__(self,other):
        return True

class Parameters:
    def __init__(self):

        self.DOBREAK = False

        # set constants

        # for timing various parts of the code
        self.last_time = 0

        # default frame rate
        self.DEFAULT_FRAME_RATE = 25.

        # for fast computation of weighted region props
        self.GRID = Grid()

        # max frac above-threshold points in a frame
        self.max_n_points_ratio = 1./250.
        # max n objects per frame to return
        self.max_n_clusters = 100 
        # for defining empty ellipses. should be obsolete, eventually
        self.empty_val = -1
        self.dummy_val = -2
        # for displaying zoom windows
        self.zoom_window_pix = 20 # 10 pixels on either side of fly center
        self.id_spinner_width = 40
        # palettes
        self.normal_palette = [[255,0,0],     # red
                               [0,255,0],     # green
                               [0,0,255],     # blue
                               [255,0,255],   # magenta
                               [255,255,0],   # yellow
                               [0,255,255],   # cyan
                               [255,127,127], # light red
                               [127,255,127], # light green
                               [127,127,255], # light blue
                               [255,127,255], # light magenta
                               [255,255,127], # light yellow
                               [127,255,255]] # light cyan
        # colorblind-friendly palette from
        # http://jfly.iam.u-tokyo.ac.jp/color/index.html
        self.colorblind_palette = [[230,159,0],   # orange
                                   [86,180,233],  # sky blue
                                   [0,158,115],   # blue-green
                                   [240,228,66],  # yellow
                                   [0,114,178],   # blue
                                   [213,94,0],    # vermillion
                                   [204,121,167]] # red-purple

        # number of frames of annotation to buffer
        self.anndata_nbuffer = 1000
        # we will store the location of every anndata_lookupinterval(th)
        # in the annotation data
        self.anndata_lookupinterval = 1000

        # background constants
        self.SHOW_BACKGROUND = 0
        self.SHOW_DISTANCE = 1
        self.SHOW_THRESH = 2
        self.SHOW_NONFORE = 3
        self.SHOW_DEV = 4
        self.SHOW_CC = 5
        self.SHOW_ELLIPSES = 6
        self.BG_TYPE_LIGHTONDARK = 0
        self.BG_TYPE_DARKONLIGHT = 1
        self.BG_TYPE_OTHER = 2
        self.BG_NORM_STD = 0
        self.BG_NORM_INTENSITY = 1
        self.BG_NORM_HOMOMORPHIC = 2
        self.ALGORITHM_MEDIAN = 0
        self.ALGORITHM_MEAN = 1
        self.print_crap = False # print debugging info
        self.watch_threads = True # print processing info
        self.watch_locks = False # print lock acquisition/releasing
        self.count_time = True # determine how long each thread works
        if not DEBUG: # turn all debugging info off
            self.print_crap = False
            self.watch_threads = False
            self.watch_locks = False
            self.count_time = False
        # display stuff
        self.status_green = "#66FF66"
        self.status_blue = "#AAAAFF"
        self.status_red = "#FF6666"
        self.status_yellow = "#FFFF66"
        self.wxvt_bg = '#DDFFDD'

        self.MAXDSHOWINFO = 10
        self.DRAW_MOTION_SCALE = 10.
        
        # set parameters to default values

        self.start_frame = 0

        # Number of identities assigned so far
        self.nids = 0

        # whether Ctrax was started in interactive mode or not
        self.interactive = True

        self.version = 0

        # number of targets
        #self.min_ntargets = 0
        #self.max_ntargets = 100

        # Movie Parameters

        # number of frames in the movie
        self.n_frames = 0
        # size of a frame
        self.movie_size = (0,0)
        # number of pixels in a frame
        self.npixels = 0
        # movie object
        self.movie = None
        self.movie_name = ''
        self.annotation_movie_name = ''

        # Background Estimation Parameters

        # number of frames to use to estimate background
        self.n_bg_frames = 100
        # which algorithm to use to estimate background
        self.use_median = True
        # interval to use for estimating background
        self.bg_firstframe = 0
        self.bg_lastframe = 99999

        # Background Subtraction Parameters

        # do we assume dark on light, light on dark, or no assumption?
        self.bg_type = self.BG_TYPE_OTHER
        # how do we normalize the background subtraction image?
        self.bg_norm_type = self.BG_NORM_STD
        # homomorphic filtering constants
        # these defaults are based on a painstaking and probably
        # inaccurate parameterization, but they're a good starting point
        self.hm_cutoff = 0.35
        self.hm_boost = 2
        self.hm_order = 2
        # number of standard deviations to threshold background at
        self.n_bg_std_thresh = 20.
        # do hysteresis
        self.n_bg_std_thresh_low = 10.
        # minimum number of standard deviations for background
        self.bg_std_min = 1.
        self.bg_std_max = 10.
        # if background intensity is greater than min_nonarena, then it
        # is not an area flies can be in, so don't allow foreground in
        # these areas
        self.min_nonarena = 256.
        # if background intensity is less than max_nonarena, then it
        # is not an area flies can be in, so don't allow foreground in
        # these areas
        self.max_nonarena = -1.
        
        # regions of interest
        self.roipolygons = []
        
        # location of arena
        self.arena_center_x = None
        self.arena_center_y = None
        self.arena_radius = None
        self.arena_edgethresh = None
        self.do_set_circular_arena = True
        # batch processing & auto-detecting arena
        self.batch_autodetect_arena = True
        self.batch_autodetect_shape = True
        self.batch_autodetect_bg_model = True
        # morphology
        self.do_use_morphology = False
        self.opening_radius = 0
        self.closing_radius = 0

        # Shape Parameters

        # upper bounds on shape
        self.maxshape = ShapeParams(10000,10000,10000,1.)
        # lower bounds on shape
        self.minshape = ShapeParams(1.,1.,1.,0.)
        self.meanshape = ShapeParams(2.64,3.56,40.25,1.98)
        self.have_computed_shape = False
        # foreground background thresh set to minbackthresh when trying to increase target area
        self.minbackthresh = 1.
        # maximum number of clusters to split a foreground connected component into during the forward pass
        self.maxclustersperblob = 5
        # maximum penalty for merging together two ccs
        self.maxpenaltymerge = 40
        # maximum area of deleted target
        self.maxareadelete = 5
        # minimum area of ignored connected components
        self.minareaignore = 2500
        # number of frames used to compute shape bounds
        self.n_frames_size = 50
        # number of standard deviations from mean for upper and lower shape bounds
        self.n_std_thresh = 3.

        # Motion model parameters

        # weight of angle in distance measure
        self.ang_dist_wt = 100.
        # maximum distance a fly can move between frames
        self.max_jump = 100.
        # dampening constant
        self.dampen = 0. # weighting term used in cvpred()
        self.angle_dampen = 0.5 

        # weight of velocity angle in choosing orientation mod 2pi
        self.velocity_angle_weight = .05
        self.max_velocity_angle_weight = .25

        # Fix errors parameters
        self.do_fix_split = True
        self.do_fix_merged = True
        self.do_fix_spurious = True
        self.do_fix_lost = True
        self.lostdetection_length = 50
        self.lostdetection_distance = 100.
        self.spuriousdetection_length = 50
        self.mergeddetection_distance = 20.
        self.mergeddetection_length = 50
        self.splitdetection_cost = 40.
        self.splitdetection_length = 50
        self.maxdcentersextra = 3.
        self.bigboundingboxextra = 2.
        # maximum number of frames to buffer for faster hindsight
        self.maxnframesbuffer = 100

        # Drawing/Display parameters

        # thickness of lines of drawn ellipses
        self.ellipse_thickness = 2
        # palette of colors to use
        self.use_colorblind_palette = False
        # colors
        self.colors = self.normal_palette
        # slider returns an integer value, so scale slider value
        # by slider resolution.
        self.sliderres = 255.

        self.tail_length = 10

        self.status_box = 0
        self.file_box = 1
        self.file_box_max_width = 40

        # how often we plot the tracking results
        self.request_refresh = False
        self.do_refresh = True
        self.framesbetweenrefresh = 1

        # Hopefully Obsolete Someday
        self.bg_est_threads = 1
        self.max_median_frames = 100
        self.max_n_obj = 500 # 50 flies can each get lost 10 times before aborting

        # Obsolete Parameters
        #self.min_cluster_size = 10 # pixels
        #self.default_thresh = 70
        #self.n_hist_bins = 100
        #self.uhistogram_cut = 5 # n-th percentile
        #self.n_frames_trunc = 2

    def __print__(self):
        s = ""
        for i,j in self.__dict__.iteritems():
            if j is None:
                s += i + ": None\n"
            else:
                s += i + ": " + str(j) + "\n"
        return s

    def __repr__( self ): return self.__print__()

    def __str__( self ): return self.__print__()

    def copy(self):
        v = Parameters()
        for i,j in self.__dict__.iteritems():
            if i == 'movie':
                continue
            try:
                v.__dict__[i] = copy.deepcopy(j)
            except:
                v.__dict__[i] = j
        return v

    def __eq__( self, other ):
        for i,j in self.__dict__.iteritems():
            if not hasattr(other,i):
                return False
            if not (j == other.__dict__[i]):
                return False
        for i,j in other.__dict__.iteritems():
            if not hasattr(self,i):
                return False
            if not (j == self.__dict__[i]):
                return False
        return True
        
params = Parameters()

class GUIConstants:
    def __init__( self ):
        self.info = AboutDialogInfo()
        self.info.SetName( "Ctrax" )
        self.info.SetVersion( __version__ )
        self.info.SetCopyright( "2007-2010, Caltech ethomics project" )
        self.info.SetDescription( """The Caltech Multiple Fly Tracker.

http://www.dickinson.caltech.edu/ctrax

Distributed under the GNU General Public License
(http://www.gnu.org/licenses/gpl.html) with
ABSOLUTELY NO WARRANTY.

This project is supported by grant R01 DA022777-01 from
the National Institute on Drug Abuse at the US NIH.""" )
        
        self.TRACK_START = "Start Tracking"
        self.TRACK_STOP = "Stop Tracking"
        self.TRACK_PLAY = "Start Playback"
        self.TRACK_PAUSE = "Stop Playback"

        params.version = __version__
# rather than global variables for each of these...

const = GUIConstants()
