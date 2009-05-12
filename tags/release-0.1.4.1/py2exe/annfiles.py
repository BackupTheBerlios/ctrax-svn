# annfiles.py
# KMB 11/06/2008

import numpy as num
import scipy.io
import wx

from version import DEBUG
from params import params
import ellipsesk as ell
import pickle

class InvalidFileFormatException(Exception):
    pass

class AnnotationFile:
    def __init__( self, filename, version, writeflag=False ):
        """This function should be called inside an exception handler
        in case the annotation header versions are different from expected
        or the file doesn't exist."""
        self.filename = filename
        self.version = version

        self.open_for_writing = False
        if writeflag:
            self.file = open( self.filename, mode="wb" )
            self.open_for_writing = True
        else:
            # check to see if file is readable and has a valid header
            self.file = open( self.filename, mode="rb" )
            self.CheckAnnHeader()
            self.file.close()

    def __del__( self ):
        if self.open_for_writing: self.file.close()

    def write_ellipses( self, ellipse_list ):
        """Write one frame of data to already-open file."""
        for ellipse in ellipse_list.itervalues():
            self.file.write( '%f\t%f\t%f\t%f\t%f\t%d\t'%(ellipse.center.x,
                                                         ellipse.center.y,
                                                         ellipse.size.width,
                                                         ellipse.size.height,
                                                         ellipse.angle,
                                                         ellipse.identity) )
        self.file.write( "\n" )

    def WriteAnnHeader( self, start_frame, background ):
        """Write the header for an annotation file."""
        SIZEOFDOUBLE = 8
        self.file.write("Ctrax header\n" )
        self.file.write("version:%s\n" %self.version )

        # parameters

        # background parameters
        self.file.write("bg_type:%d\n"%params.bg_type)
        self.file.write("n_bg_std_thresh:%.1f\n" %params.n_bg_std_thresh )
        self.file.write("n_bg_std_thresh_low:%.1f\n" %params.n_bg_std_thresh_low )
        self.file.write("bg_std_min:%.2f\n" %params.bg_std_min)
        self.file.write("bg_std_max:%.2f\n" %params.bg_std_max)
        self.file.write("n_bg_frames:%d\n" %params.n_bg_frames)
        self.file.write("min_nonarena:%.1f\n" %params.min_nonarena)
        self.file.write("max_nonarena:%.1f\n" %params.max_nonarena)
        if params.arena_center_x is not None:
            self.file.write("arena_center_x:%.2f\n"%params.arena_center_x)
            self.file.write("arena_center_y:%.2f\n"%params.arena_center_y)
            self.file.write("arena_radius:%.2f\n"%params.arena_radius)
        self.file.write("do_set_circular_arena:%d\n"%params.do_set_circular_arena)
        self.file.write("do_use_morphology:%d\n" %params.do_use_morphology)
        self.file.write("opening_radius:%d\n" %params.opening_radius)
        self.file.write("closing_radius:%d\n" %params.closing_radius)
        self.file.write("bg_algorithm:")
        if params.use_median:
            self.file.write("median\n")
        else:
            self.file.write("mean\n")
        if hasattr(background,'med'):
            sz = background.med.size*SIZEOFDOUBLE
            self.file.write("background median:%d\n"%sz)
            self.file.write(background.med)
        if hasattr(background,'mean'):
            sz = background.mean.size*SIZEOFDOUBLE
            self.file.write("background mean:%d\n"%sz)
            self.file.write(background.mean)
        self.file.write('bg_norm_type:%d\n'%params.bg_norm_type)
        if hasattr(background,'mad'):
            sz = background.mad.size*SIZEOFDOUBLE
            self.file.write("background mad:%d\n"%sz)
            self.file.write(background.mad)
        if hasattr(background,'std'):
            sz = background.std.size*SIZEOFDOUBLE
            self.file.write("background std:%d\n"%sz)
            self.file.write(background.std)
        if hasattr(background,'hfnorm'):
            sz = background.hfnorm.size*SIZEOFDOUBLE
            self.file.write("hfnorm:%d\n"%sz)
            self.file.write(background.hfnorm)
        if hasattr(params,'roipolygons'):
            s = pickle.dumps(params.roipolygons)
            self.file.write("roipolygons:%d\n"%len(s))
            self.file.write(s)
        self.file.write('bg_norm_type:%d\n'%params.bg_norm_type)
        self.file.write('hm_cutoff:%.2f\n'%params.hm_cutoff)
        self.file.write('hm_boost:%d\n'%params.hm_boost)
        self.file.write('hm_order:%d\n'%params.hm_order)

        # shape parameters
        self.file.write('maxarea:%.2f\n'%params.maxshape.area )
        self.file.write('maxmajor:%.2f\n'%params.maxshape.major )
        self.file.write('maxminor:%.2f\n'%params.maxshape.minor )
        self.file.write('maxecc:%.2f\n'%params.maxshape.ecc )
        self.file.write('minarea:%.2f\n'%params.minshape.area )
        self.file.write('minmajor:%.2f\n'%params.minshape.major )
        self.file.write('minminor:%.2f\n'%params.minshape.minor )
        self.file.write('minecc:%.2f\n'%params.minshape.ecc )
        self.file.write('meanarea:%.2f\n'%params.meanshape.area )
        self.file.write('meanmajor:%.2f\n'%params.meanshape.major )
        self.file.write('meanminor:%.2f\n'%params.meanshape.minor )
        self.file.write('meanecc:%.2f\n'%params.meanshape.ecc )
        self.file.write('nframes_size:%d\n'%params.n_frames_size)
        self.file.write('nstd_shape:%d\n'%params.n_std_thresh)

        # motion model parameters
        self.file.write('max_jump:%.2f\n'%params.max_jump)
        self.file.write('ang_dist_wt:%.2f\n'%params.ang_dist_wt)
        self.file.write('center_dampen:%.2f\n'%params.dampen)
        self.file.write('angle_dampen:%.2f\n'%params.angle_dampen)
        if params.velocity_angle_weight is not None:
            self.file.write('velocity_angle_weight:%.2f\n'%params.velocity_angle_weight)
        if params.max_velocity_angle_weight is not None:
            self.file.write('max_velocity_angle_weight:%.2f\n'%params.max_velocity_angle_weight)

        # observation parameters
        self.file.write('minbackthresh:%.2f\n'%params.minbackthresh)
        self.file.write('maxpenaltymerge:%.2f\n'%params.maxpenaltymerge)
        self.file.write('maxareadelete:%.2f\n'%params.maxareadelete)
        self.file.write('minareaignore:%.2f\n'%params.minareaignore)

        # hindsight parameters
        self.file.write('do_fix_split:%d\n'%params.do_fix_split)
        self.file.write('splitdetection_length:%d\n'%params.splitdetection_length)
        self.file.write('splitdetection_cost:%.2f\n'%params.splitdetection_cost)

        self.file.write('do_fix_merged:%d\n'%params.do_fix_merged)
        self.file.write('mergeddetection_length:%d\n'%params.mergeddetection_length)
        self.file.write('mergeddetection_distance:%.2f\n'%params.mergeddetection_distance)
        
        self.file.write('do_fix_spurious:%d\n'%params.do_fix_spurious)
        self.file.write('spuriousdetection_length:%d\n'%params.spuriousdetection_length)

        self.file.write('do_fix_lost:%d\n'%params.do_fix_lost)
        self.file.write('lostdetection_length:%d\n'%params.lostdetection_length)

        self.file.write('movie_name:' + params.movie_name + '\n')

        self.max_jump = params.max_jump
        self.file.write('start_frame:%d\n'%start_frame)
        
        self.file.write('data format:identity x y major minor angle\n')

        self.file.write('end header\n')

        self.file.flush()

        self.endofheader = self.file.tell()

    def CheckAnnHeader( self ):
        """Read header from an annotation file."""

        self.file.seek(0,0)

        # first line says "Ctrax header\n" or "mtrax header\n"
        # note that Ctrax used to be called mtrax, and we have left 
        # the annotation file format the same
        line = self.file.readline()
        if line == 'mtrax header\n':
            pass
        elif line == 'Ctrax header\n':
            pass
        else:
            print "line = >" + line + "<"
            raise InvalidFileFormatException("Annotation file does not start with 'Ctrax header' or 'mtrax header'")
        
        # next lines have parameters.
        i = 0
        while True:
            line = self.file.readline()
            # header ends when we see 'end header\n'
            if line == '':
                raise InvalidFileFormatException("End of Annotation File reached; did not find 'end header'")
            if line == 'end header\n':
                self.start_data = self.file.tell()
                break
            # split into parameter and value at :
            words = line.split(':',1)
            if len(words) is not 2:
                raise InvalidFileFormatException("More than one ':' in line >> " + line)
            parameter = words[0]
            value = words[1]
            # strip the \n
            if value[-1] == '\n':
                value = value[:-1]
            else:
                raise InvalidFileFormatException("Line does not end in newline character. line >> " + line)
            if len(value) == 0:
                continue
            if parameter == 'background median' or parameter == 'background mean' or \
                   parameter == 'background mad' or parameter == 'background std' or \
                   parameter == 'hfnorm' or parameter == 'roipolygons':
                sz = int(value)
                self.file.seek(sz,1)
            elif parameter == 'data format':
                self.n_fields = len(value.split())

            i += 1

    def ReadAnnHeader( self, background ):
        """Read header from an annotation file."""

        # first line says "mtrax header\n"
        line = self.file.readline()
        if line == 'mtrax header\n':
            pass
        elif line == 'Ctrax header\n':
            pass
        else:
            print "line = >" + line + "<"
            raise InvalidFileFormatException("Annotation file does not start with 'Ctrax header' or 'mtrax header'")
        
        # next lines have parameters.
        while True:
            line = self.file.readline()
            # header ends when we see 'end header\n'
            if line == '':
                raise InvalidFileFormatException("End of Annotation File reached; did not find 'end header'")
            if line == 'end header\n':
                self.start_data = self.file.tell()
                break
            # split into parameter and value at :
            words = line.split(':',1)
            if len(words) is not 2:
                raise InvalidFileFormatException("More than one ':' in line >> " + line)
            parameter = words[0]
            value = words[1]
            # strip the \n
            if value[-1] == '\n':
                value = value[:-1]
            else:
                raise InvalidFileFormatException("Line does not end in newline character. line >> " + line)
            if len(value) == 0:
                continue
            if parameter == 'bg_type':
                params.bg_type = float(value)
            elif parameter == 'n_bg_std_thresh':
                params.n_bg_std_thresh = float(value)
            elif parameter == 'n_bg_std_thresh_low':
                params.n_bg_std_thresh_low = float(value)
            elif parameter == 'bg_std_min':
                params.bg_std_min = float(value)
            elif parameter == 'bg_std_max':
                params.bg_std_max = float(value)
            elif parameter == 'n_bg_frames':
                params.n_bg_frames = int(value)
            elif parameter == 'min_nonarena':
                params.min_nonarena = float(value)                
            elif parameter == 'max_nonarena':
                params.max_nonarena = float(value)
            elif parameter == 'arena_center_x':
                params.arena_center_x = float(value)
            elif parameter == 'arena_center_y':
                params.arena_center_y = float(value)
            elif parameter == 'arena_radius':
                params.arena_radius = float(value)
            elif parameter == 'do_set_circular_arena':
                params.do_set_circular_arena = num.bool(value)
            elif parameter == 'do_use_morphology':
                params.do_use_morphology = num.bool(value)
            elif parameter == 'opening_radius':
                params.opening_radius = int(value)
                background.opening_struct = background.create_morph_struct(params.opening_radius)
            elif parameter == 'closing_radius':
                params.closing_radius = int(value)
                background.closing_struct = background.create_morph_struct(params.closing_radius)
            elif parameter == 'bg_algorithm':
                if value == 'median':
                    params.use_median = True
                else:
                    params.use_median = False
            elif parameter == 'background median':
                sz = int(value)
                background.med = num.fromstring(self.file.read(sz),'<d')
                background.med.shape = params.movie_size
            elif parameter == 'background mean':
                sz = int(value)
                background.mean = num.fromstring(self.file.read(sz),'<d')
                background.mean.shape = params.movie_size
            elif parameter == 'bg_norm_type':
                params.bg_norm_type = int(value)
            elif parameter == 'background mad':
                sz = int(value)
                background.mad = num.fromstring(self.file.read(sz),'<d')
                background.mad.shape = params.movie_size
            elif parameter == 'background std':
                sz = int(value)
                background.std = num.fromstring(self.file.read(sz),'<d')
                background.std.shape = params.movie_size
            elif parameter == 'hfnorm':
                sz = int(value)
                background.hfnorm = num.fromstring(self.file.read(sz),'<d')
                background.hfnorm.shape = params.movie_size
            elif parameter == 'roipolygons':
                sz = int(value)
                params.roipolygons = pickle.loads(self.file.read(sz))
            elif parameter == 'hm_cutoff':
                params.hm_cutoff = float(value)
            elif parameter == 'hm_boost':
                params.hm_boost = int(value)
            elif parameter == 'hm_order':
                params.hm_order = int(value)
            elif parameter == 'maxarea':
                params.maxshape.area = float(value)
            elif parameter == 'maxmajor':
                params.maxshape.major = float(value)
            elif parameter == 'maxminor':
                params.maxshape.minor = float(value)
            elif parameter == 'maxecc':
                params.maxshape.ecc = float(value)
            elif parameter == 'minarea':
                params.minshape.area = float(value)
            elif parameter == 'minmajor':
                params.minshape.major = float(value)
            elif parameter == 'minminor':
                params.minshape.minor = float(value)
            elif parameter == 'minecc':
                params.minshape.ecc = float(value)
            elif parameter == 'meanarea':
                params.meanshape.area = float(value)
            elif parameter == 'meanmajor':
                params.meanshape.major = float(value)
            elif parameter == 'meanminor':
                params.meanshape.minor = float(value)
            elif parameter == 'meanecc':
                params.meanshape.ecc = float(value)
            elif parameter == 'nframes_size':
                params.n_frames_size = int(value)
            elif parameter == 'nstd_shape':
                params.n_std_thresh = float(value)
            elif parameter == 'max_jump':
                params.max_jump = float(value)
                self.max_jump = params.max_jump
            elif parameter == 'ang_dist_wt':
                params.ang_dist_wt = float(value)
            elif parameter == 'center_dampen':
                params.dampen = float(value)
            elif parameter == 'angle_dampen':
                params.angle_dampen = float(value)
            elif parameter == 'velocity_angle_weight':
                params.velocity_angle_weight = float(value)
            elif parameter == 'max_velocity_angle_weight':
                params.max_velocity_angle_weight = float(value)
            elif parameter == 'minbackthresh':
                params.minbackthresh = float(value)
            elif parameter == 'maxpenaltymerge':
                params.maxpenaltymerge = float(value)
            elif parameter == 'maxareadelete':
                params.maxareadelete = float(value)
            elif parameter == 'minareaignore':
                params.minareaignore = float(value)
            elif parameter == 'do_fix_split':
                params.do_fix_split = bool(value)
            elif parameter == 'splitdetection_length':
                params.splitdetection_length = int(value)
            elif parameter == 'splitdetection_cost':
                params.splitdetection_cost = float(value)
            elif parameter == 'do_fix_merged':
                params.do_fix_merged = bool(value)
            elif parameter == 'mergeddetection_length':
                params.mergeddetection_length = int(value)
            elif parameter == 'mergeddetection_distance':
                params.mergeddetection_distance = float(value)
            elif parameter == 'do_fix_spurious':
                params.do_fix_spurious = bool(value)
            elif parameter == 'spuriousdetection_length':
                params.spuriousdetection_length = int(value)
            elif parameter == 'do_fix_lost':
                params.do_fix_lost = bool(value)
            elif parameter == 'lostdetection_length':
                params.lostdetection_length = int(value)
            elif parameter == 'movie_name':
                params.annotation_movie_name = value
            elif parameter == 'start_frame':
                params.start_frame = int(value)
            elif parameter == 'data format':
                self.n_fields = len(value.split())

        background.updateParameters()

    def ReadSettings( self ):
        """Read header from an annotation file."""

        # first line says "mtrax header\n"
        line = self.file.readline() 
        if line == 'mtrax header\n':
            pass
        elif line == 'Ctrax header\n':
            pass
        else:
            print "line = >" + line + "<"
            raise InvalidFileFormatException("Annotation file does not start with 'Ctrax header' or 'mtrax header'")
        
        # next lines have parameters.
        while True:
            line = self.file.readline()
            # header ends when we see 'end header\n'
            if line == '':
                raise InvalidFileFormatException("End of Annotation File reached; did not find 'end header'")
            if line == 'end header\n':
                self.start_data = self.file.tell()
                break
            # split into parameter and value at :
            words = line.split(':',1)
            if len(words) is not 2:
                raise InvalidFileFormatException("More than one ':' in line >> " + line)
            parameter = words[0]
            value = words[1]
            # strip the \n
            if value[-1] == '\n':
                value = value[:-1]
            else:
                raise InvalidFileFormatException("Line does not end in newline character. line >> " + line)
            if len(value) == 0:
                continue
            if parameter == 'bg_type':
                params.bg_type = float(value)
            elif parameter == 'n_bg_std_thresh':
                params.n_bg_std_thresh = float(value)
            elif parameter == 'n_bg_std_thresh_low':
                params.n_bg_std_thresh_low = float(value)
            elif parameter == 'bg_std_min':
                params.bg_std_min = float(value)
            elif parameter == 'bg_std_max':
                params.bg_std_max = float(value)
            elif parameter == 'n_bg_frames':
                params.n_bg_frames = int(value)
            elif parameter == 'min_nonarena':
                params.min_nonarena = float(value)                
            elif parameter == 'max_nonarena':
                params.max_nonarena = float(value)
            elif parameter == 'arena_center_x':
                params.arena_center_x = float(value)
            elif parameter == 'arena_center_y':
                params.arena_center_y = float(value)
            elif parameter == 'arena_radius':
                params.arena_radius = float(value)
            elif parameter == 'do_set_circular_arena':
                params.do_set_circular_arena = num.bool(value)
            elif parameter == 'do_use_morphology':
                params.do_use_morphology = num.bool(value)
            elif parameter == 'opening_radius':
                params.opening_radius = int(value)
            elif parameter == 'closing_radius':
                params.closing_radius = int(value)
            elif parameter == 'bg_algorithm':
                if value == 'median':
                    params.use_median = True
                else:
                    params.use_median = False
            # note that when loading settings, the background model is not
            # read in
            elif parameter == 'background median':
                sz = int(value)
                self.file.seek(sz,1)
            elif parameter == 'background mean':
                sz = int(value)
                self.file.seek(sz,1)
            elif parameter == 'bg_norm_type':
                params.bg_norm_type = int(value)
            elif parameter == 'background mad':
                sz = int(value)
                self.file.seek(sz,1)
            elif parameter == 'background std':
                sz = int(value)
                self.file.seek(sz,1)
            elif parameter == 'hfnorm':
                sz = int(value)
                self.file.seek(sz,1)
            elif parameter == 'roipolygons':
                sz = int(value)
                self.file.seek(sz,1)
            elif parameter == 'hm_cutoff':
                params.hm_cutoff = float(value)
            elif parameter == 'hm_boost':
                params.hm_boost = int(value)
            elif parameter == 'hm_order':
                params.hm_order = int(value)
            elif parameter == 'maxarea':
                params.maxshape.area = float(value)
            elif parameter == 'maxmajor':
                params.maxshape.major = float(value)
            elif parameter == 'maxminor':
                params.maxshape.minor = float(value)
            elif parameter == 'maxecc':
                params.maxshape.ecc = float(value)
            elif parameter == 'minarea':
                params.minshape.area = float(value)
            elif parameter == 'minmajor':
                params.minshape.major = float(value)
            elif parameter == 'minminor':
                params.minshape.minor = float(value)
            elif parameter == 'minecc':
                params.minshape.ecc = float(value)
            elif parameter == 'meanarea':
                params.meanshape.area = float(value)
            elif parameter == 'meanmajor':
                params.meanshape.major = float(value)
            elif parameter == 'meanminor':
                params.meanshape.minor = float(value)
            elif parameter == 'meanecc':
                params.meanshape.ecc = float(value)
            elif parameter == 'nframes_size':
                params.n_frames_size = int(value)
            elif parameter == 'nstd_shape':
                params.n_std_thresh = float(value)
            elif parameter == 'max_jump':
                params.max_jump = float(value)
                self.max_jump = params.max_jump
            elif parameter == 'ang_dist_wt':
                params.ang_dist_wt = float(value)
            elif parameter == 'center_dampen':
                params.dampen = float(value)
            elif parameter == 'angle_dampen':
                params.angle_dampen = float(value)
            elif parameter == 'velocity_angle_weight':
                params.velocity_angle_weight = float(value)
            elif parameter == 'max_velocity_angle_weight':
                params.max_velocity_angle_weight = float(value)
            elif parameter == 'minbackthresh':
                params.minbackthresh = float(value)
            elif parameter == 'maxpenaltymerge':
                params.maxpenaltymerge = float(value)
            elif parameter == 'maxareadelete':
                params.maxareadelete = float(value)
            elif parameter == 'minareaignore':
                params.minareaignore = float(value)
            elif parameter == 'data format':
                self.n_fields = len(value.split())
            elif parameter == 'movie_name':
                params.annotation_movie_name = value

    def ParseData( self, line):
        """Split a line of annotation data into per-fly information.
        Returns an EllipseList. Allows passing in an EllipseList, which is
        overwritten and returned, to avoid memory reallocation."""
        
        fly_sep = line.split()
        
        # initialize
        ellipses = ell.TargetList()

        for ff in range( len(fly_sep)/self.n_fields ):
            # parse data from line
            new_ellipse = ell.Ellipse( centerX=float(fly_sep[self.n_fields*ff]),
                                       centerY=float(fly_sep[self.n_fields*ff+1]),
                                       sizeW=float(fly_sep[self.n_fields*ff+2]),
                                       sizeH=float(fly_sep[self.n_fields*ff+3]),
                                       angle=float(fly_sep[self.n_fields*ff+4]),
                                       identity=int(fly_sep[self.n_fields*ff+5]) )
            if new_ellipse.identity >= params.nids:
                params.nids = new_ellipse.identity+1

            ellipses.append( new_ellipse )

        return ellipses
    
    def GetAnnotation( self, interactive, background ):
        """Read all frames of annotation data. Returns a list of TargetLists."""
        # open file and grab header
        self.file = open( self.filename, mode="rb" )

        self.ReadAnnHeader(background)

        ellipses = []

        for i, line in enumerate( self.file ):
            ellipses.append(self.ParseData(line))
            if interactive: wx.Yield()
        self.file.close()

        return ellipses

    def WriteMAT( self, filename, data, dosavestamps=False ):
        """Writes a MATLAB .mat file from the data sent in. Expects a
        list of EllipseLists."""
        # how many targets per frame
        nframes = len(data)
        startframe = params.start_frame
        ntargets = num.ones( nframes )
        for i,frame in enumerate(data):
            ntargets[i] = len(frame)
        z = num.sum(ntargets)

        # allocate arrays as all NaNs
        x_pos = num.ones( z ) * num.nan
        y_pos = num.ones( z ) * num.nan
        maj_ax = num.ones( z ) * num.nan
        min_ax = num.ones( z ) * num.nan
        angle = num.ones( z ) * num.nan
        identity = num.ones( z ) * num.nan

        # store data in arrays
        i = 0
        for ff in range( nframes ):
            for id, ee in data[ff].iteritems():
                x_pos[i] = ee.center.x
                y_pos[i] = ee.center.y
                maj_ax[i] = ee.height
                min_ax[i] = ee.width
                angle[i] = ee.angle
                identity[i] = ee.identity
                i+=1

        # save data to mat-file
        save_dict = {'x_pos': x_pos,
                     'y_pos': y_pos,
                     'maj_ax': maj_ax,
                     'min_ax': min_ax,
                     'angle': angle,
                     'identity': identity,
                     'ntargets': ntargets,
                     'startframe': startframe}

        if dosavestamps:
            print "reading in movie to save stamps"
            stamps = params.movie.get_some_timestamps(t1=startframe,t2=startframe+nframes)
            save_dict['timestamps'] = stamps

        scipy.io.savemat( filename, save_dict )
        # menu_file_write_timestamps
    def RewriteTracks(self,data):

        # open the file for reading and writing
        self.file = open( self.filename, mode="rb+" )

        # find the end of the header
        self.CheckAnnHeader()

        # truncate here
        self.file.truncate()
        
        # write the data
        for ff in range(len(data)):
            self.write_ellipses(data[ff])

        self.file.close()

def RewriteTracks(filename,data):

    # create the annotation file structure
    annfile = AnnotationFile(filename,params.version,writeflag=False)

    # main function
    annfile.RewriteTracks(data)
