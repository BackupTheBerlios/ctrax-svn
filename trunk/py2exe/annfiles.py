# annfiles.py
# KMB 11/06/2008

import numpy as num
import scipy.io
import wx
import os
import tempfile
import shutil

from version import __version__
from version import DEBUG_ANN as DEBUG
from params import params
import ellipsesk as ell
import pickle

dataformatstring = 'identity x y major minor angle'

class InvalidFileFormatException(Exception):
    pass

class AnnotationFile:

    def __init__( self, filename=None, bg_imgs=None, doreadheader=True,
                  justreadheader=False, doreadbgmodel=True ):
        """This function should be called inside an exception handler
        in case the annotation header versions are different from expected
        or the file doesn't exist."""

        # name of annotation file
        self.filename = filename

        # version of annotation file
        self.version = __version__

        # file pointer
        self.file = None

        # background data
        self.bg_imgs = bg_imgs

        # read in the header if that is necessary (filename should exist)
        if justreadheader:
            if DEBUG: print "Just reading header"
            self.file = open(filename,'rb')
            if doreadbgmodel:
                if DEBUG: print "Reading header ..."
                self.ReadAnnHeader()
                if DEBUG: print "Done reading header."
            else:
                if DEBUG: print "Reading settings ..."
                self.ReadSettings()
                if DEBUG: print "Done reading settings."
            self.file.close()
            return

        # initialize data structures that no tracks have been read
        self.isdatawritten = False
        self.firstframetracked = 0
        self.lastframetracked = -1
        self.nframestracked = 0
        self.firstframewritten = 0
        self.lastframewritten = -1
        self.firstframebuffered = 0
        self.lastframebuffered = -1
        self.buffer = []
        self.n = 0
        self.nbuffer = params.anndata_nbuffer
        self.lookupinterval = params.anndata_lookupinterval
        self.lookup = []
        self.isdatawritten = False
        self.idtable = dict()
        self.n_fields = len(dataformatstring.split())
        self.recycledids = []

        if DEBUG: print "Initialized data structures"

        # create a temporary annotation file
        if filename is None:
            if DEBUG: print "Creating a temporary file"
            # create a temporary file
            self.file = tempfile.NamedTemporaryFile(mode='wb+',suffix='.ann')
            self.filename = self.file.name
            if DEBUG: print "temporary filename = " + self.filename
            params.nids = 0
            return

        # check if filename exists
##        if not os.path.isfile(filename):
##            # if not, then create the file
##            print "Annotation file %s does not exist yet, not opening existing annotation"%filename
##            self.file = open(filename,'wb+')
##            params.nids = 0
##            return

##        if DEBUG: print "Annotation file exists"
        # checked below instead... allow additional parameters to be set

        # open the file for reading
        newfile = False
        try:
            self.file = open(filename,'rb+')
        except IOError:
            self.file = open(filename,'wb+')
            newfile = True

        # check the annotation header
        if not newfile:
            try:
                if DEBUG: print "Checking annotation header..."
                self.CheckAnnHeader()
                if DEBUG: print "Done checking annotation header."
            except InvalidFileFormatException:
                if params.interactive and not params.batch_executing:
                    # make sure it is okay to overwrite this file
                    msgtxt = "Annotation file %s exists, but could not parse the header. Overwrite?"%filename
                    if wx.MessageBox(msgtxt,"Overwrite %s?"%filename,wx.OK|wx.CANCEL) == wx.CANCEL:
                        self.filename = None
                        # do not reset params.nids
                        return
                self.file.close()
                # overwrite
                self.file = open(filename,'wb+')
                newfile = True

        params.nids = 0

        # update nbuffer
        # maximum number of frames we'll need to consider at a time
        self.maxlookback = max(params.lostdetection_length,
                               params.spuriousdetection_length,
                               params.mergeddetection_length,
                               params.splitdetection_length)

        # number of frames to buffer
        self.nbuffer = max(params.anndata_nbuffer,self.maxlookback)
        self.lookupinterval = params.anndata_lookupinterval

        if DEBUG: print "nbuffer = %d, lookupinterval = %d"%(self.nbuffer,self.lookupinterval)

        # in non-interactive mode or for brand-new annotation files,
        # we will not read in the trajectories
        if (not params.interactive) or params.batch_executing or newfile:
            return

        # read in the header
        if doreadheader:
            self.ReadAnnHeader()
        # if we don't read it, we will still be at the correct file location from call to CheckAnnHeader
        
        # check if we want to read in the trajectories from this file
        msgtxt = "Read old trajectories from %s? Choosing not to read will cause the trajectories in this file to be overwritten."%filename
        if wx.MessageBox(msgtxt,"Read trajectories?",wx.YES|wx.NO) == wx.NO:
            if DEBUG: print "Not reading trajectories"
            return


        if DEBUG: print "Reading trajectories..."

        # first frame tracked is params.start_frame
        self.firstframetracked = params.start_frame
        self.firstframebuffered = self.firstframetracked
        # have not read anything in yet
        self.lastframetracked = self.firstframetracked-1
        self.nframestracked = 0
        self.lastframebuffered = self.firstframebuffered-1
        self.n = 0

        if DEBUG: print "initialized firstframetracked = %d, firstframebuffered = %d"%(self.firstframetracked,self.firstframebuffered)

        i = 0
        while True:

            # read line and add to lookup if necessary
            p = self.file.tell()
            line = self.file.readline()
            if line == '':
                break
            if i%self.lookupinterval == 0:
                if DEBUG: print "i = %d, adding %x to the lookup table for frame = %d"%(i,p,self.lastframetracked+1)
                if DEBUG: print "line = " + str(line.split())
                self.lookup.append(p)

            ells = self.ParseData(line)
            # count ids
            self.CountIds(ells)
            #ells = self.ReplaceIds(ells)

            #print "ells[%d] = "%i + str(ells)

            # add to the buffer if one of the first nbuffer frames
            if self.n < self.nbuffer:
                self.buffer.append(ells)
                self.lastframebuffered+=1
                self.n+=1

            # increment
            i+=1
            self.lastframetracked+=1

        self.nframestracked = self.lastframetracked - self.firstframetracked + 1

        if DEBUG: print "Finished reading in trajectories"
        if DEBUG: print "firstframetracked = %d, lastframetracked = %d"%(self.firstframetracked,self.lastframetracked)
        if DEBUG: print "firstframebuffered = %d, lastframebuffered = %d, n = %d"%(self.firstframebuffered,self.lastframebuffered,self.n)
        # frames written = frames tracked
        self.firstframewritten = self.firstframetracked
        self.lastframewritten = self.lastframetracked

        self.isdatawritten = True

        if DEBUG: print "Finished reading trajectories"
        if DEBUG: print "framestracked = [%d,%d], framesbuffered = [%d,%d], frameswritten = [%d,%d]"%(self.firstframetracked,self.lastframetracked,self.firstframebuffered,self.lastframebuffered,self.firstframewritten,self.lastframewritten)

        return

    def InitializeData(self,firstframe=0,lastframe=None):
        """Initialize data structures, annotation file to store 
        firstframe through lastframe.
        Data structures:
        buffer: a queue that will store up to maxnbuffer frames worth of 
        annotation. It is a list of lists of ellipses. the first frame
        stored is "firstframebuffered" and the last frame stored is 
        "lastframebuffered". the number of frames stored so far is "n" and
        at most "nbuffer" will be stored. 
        file: the annotation file contains the header followed by the 
        ellipses for each frame from "firstframewritten" through
        "lastframewritten". 
        lookup: list of the locations in the file of tracks for
        certain frames. we store the locations of frames f such that
        (f - firstframetracked) % lookupinterval == 0
        """

        if lastframe is None:
            lastframe = firstframe - 1

        # make sure that if lastframe is >= firstframe, data is written
        if (lastframe >= firstframe) and (not self.isdatawritten):
            raise Exception, 'This should never happen: we are trying to store frames, but we have not written the annotation file yet'
        
        if (lastframe < firstframe):
            self.InitializeEmptyData(firstframe,lastframe)
        else:
            self.CropData(firstframe,lastframe)

        if DEBUG: print "Finished Initializing Data"

    def InitializeEmptyData(self,firstframe,lastframe):

        self.firstframetracked = firstframe
        self.lastframetracked = lastframe
        self.nframestracked = 0
        self.firstframewritten = firstframe
        self.lastframewritten = lastframe
        self.firstframebuffered = firstframe
        self.lastframebuffered = lastframe
        self.buffer = []
        self.n = 0
        self.lookup = []
        self.lookupinterval = params.anndata_lookupinterval
        self.idtable = dict()
        self.recycledids = []
        params.nids = 0
        
        # maximum number of frames we'll need to consider at a time
        self.maxlookback = max(params.lostdetection_length,
                               params.spuriousdetection_length,
                               params.mergeddetection_length,
                               params.splitdetection_length)

        # number of frames to buffer
        self.nbuffer = max(params.anndata_nbuffer,self.maxlookback)

        if not self.file.closed:
            self.file.close()

        # open the file and write the annotation header
        self.file = open( self.filename, mode="wb+" )

        # write the header
        if DEBUG: print "Writing annotation header"
        self.WriteAnnHeader(params.start_frame)

    def InitializeBufferForTracking(self,f):
        if not self.IsAnnData() or f <= self.firstframetracked:
            return
        
        g = max(f-self.nbuffer,self.firstframetracked)
        if DEBUG: print "Filling buffer with frames %d through %d"%(g,f-1)
        self.seek(g)
        self.buffer = []
        self.firstframebuffered = g
        self.lastframebuffered = f-1
        self.n = self.lastframebuffered - self.firstframebuffered + 1
        for h in range(g,f):
            self.buffer.append(self.read_ellipses())

        # we will be rewriting the previous maxlookback frames
        g = max(self.firstframetracked,f-self.maxlookback)
        self.truncate(g)

    def truncate(self,g):

        # truncate the ann file
        self.seek(g)
        if DEBUG: print "Cropping written trajectories before frame %d = %x"%(g,self.file.tell())
        self.file.truncate()

        # may also need to truncate lookup: i is the last lookup that 
        # we should keep
        (i,j) = self.lookupfloor(g-1)
        for j in range(i+1,len(self.lookup)):
            if DEBUG: print "Removing lookup[%d] = %x from lookup"%(j,self.lookup[j])
            tmp = self.lookup.pop(j)

        # update lastframewritten
        self.lastframewritten = g-1


    def CropData(self,firstframe,lastframe):

        # sanity check: we should only be copying frames that we 
        # have stored in the annotation file
        if firstframe < self.firstframewritten or \
                lastframe > self.lastframewritten:
            raise Exception, 'This should never happen: we need to store frames that we do not have in the annotation file'

        if self.file.closed:
            self.file = open( self.filename, mode="rb+" )
        else:
            # otherwise, make sure all writing is done and 
            # seek to the start
            self.file.flush()

        # create a temporary file and copy old file to this 
        # temporary file
        self.file.seek(0,0)
        cpfile = tempfile.TemporaryFile()
        shutil.copyfileobj(self.file,cpfile)
        cpfile.flush()

        if DEBUG: print "Copied ann file to temporary file"

        # close the original file
        self.file.close()
        
        # open file for writing
        self.file = open( self.filename, mode="wb+" )

        if DEBUG: print "Reopened %s for writing"%self.filename

        # write the header
        self.WriteAnnHeader(params.start_frame)

        if DEBUG: print "Wrote header"

        # copy the old lookup data structures
        oldlookup = self.lookup
        oldlookupinterval = self.lookupinterval
        oldfirstframetracked = self.firstframetracked

        # initialize the new lookup data structures
        
        # lookup table currently empty
        self.lookup = []
        # grab lookupinterval from parameters
        self.lookupinterval = params.anndata_lookupinterval

        # start tracking and writing at firstframe
        self.firstframetracked = firstframe
        self.firstframewritten = firstframe

        # maximum number of frames we'll need to consider at a time
        self.maxlookback = max(params.lostdetection_length,
                               params.spuriousdetection_length,
                               params.mergeddetection_length,
                               params.splitdetection_length)

        # stop tracking so far at lastframe
        self.lastframetracked = lastframe
        # stop writing at lastframe
        self.lastframewritten = lastframe

        # number of frames to buffer
        self.nbuffer = max(params.anndata_nbuffer,self.maxlookback)

        # we will buffer the initial frames, since we will 
        # create hindsight structure, then call InitializeBuffer
        # before in-order tracking
        self.firstframebuffered = firstframe
        self.lastframebuffered = min(lastframe,firstframe+self.nbuffer-1)
        self.n = self.lastframebuffered - self.firstframebuffered + 1
        self.buffer = []

        if DEBUG: print "CropData: firstframe = %d, lastframe = %d, framestracked = [%d,%d], framesbuffered = [%d,%d], frameswritten = [%d,%d]"%(firstframe,lastframe,self.firstframetracked,self.lastframetracked,self.firstframebuffered,self.lastframebuffered,self.firstframewritten,self.lastframewritten)

        # seek to frame firstframe in the old annotation file
        self.seek(firstframe,
                  filep=cpfile,
                  lookup=oldlookup,
                  lookupinterval=oldlookupinterval,
                  firstframetracked=oldfirstframetracked)
            
        # go through all the frames to copy
        params.nids = 0
        self.idtable = {}
        for f in range(firstframe,lastframe+1):

            # add to new lookup table if necessary
            if self.islookupframe(f):
                self.file.flush()
                self.lookup.append(self.file.tell())
                if DEBUG: print "Adding frame %d (%x) to new lookup"%(f,self.lookup[-1])

            ells = self.read_ellipses(filein=cpfile)
            # count and REPLACE ids
            ells = self.ReplaceIds(ells)

            if DEBUG and self.islookupframe(f):
                print "writing frame %d at %x: "%(f,self.file.tell()) + str(ells)

            # add to buffer if necessary
            if (f >= self.firstframebuffered) and \
                    (f <= self.lastframebuffered):
                self.buffer.append(ells)

            # copy to new file
            self.write_ellipses(ells)
            #self.file.write(s)
            
        # close (and delete automatically) copy
        cpfile.close()

        self.nframestracked = self.lastframetracked - self.firstframetracked + 1

        self.isdatawritten = True
        if DEBUG: print "Finished cropping"

    def seek(self,f,filep=None,lookup=None,lookupinterval=None,firstframetracked=None):

        if filep is None:
            filep = self.file
        if lookup is None:
            lookup = self.lookup
        if lookupinterval is None:
            lookupinterval = self.lookupinterval
        if firstframetracked is None:
            firstframetracked = self.firstframetracked
        
        (i,g) = self.lookupfloor(f,lookupinterval,
                                 firstframetracked)
        if DEBUG: print "seek: lookup[%d] = %x + %d lines"%(i,lookup[i],f-g)
        filep.seek(lookup[i],0)
        for h in range(g,f):
            tmp = filep.readline()
        if DEBUG: print "Frame %d at %x"%(f,filep.tell())

        
    def islookupframe(self,f,lookupinterval=None,firstframetracked=None):
        if lookupinterval is None:
            lookupinterval = self.lookupinterval
        if firstframetracked is None:
            firstframetracked = self.firstframetracked

        return 0 == ((f - firstframetracked) % lookupinterval)

    def get_frame(self,f):

        # refill buffer if necessary
        #if DEBUG: print "get_frame %d, framestracked = [%d,%d], framesbuffered = [%d,%d], frameswritten = [%d,%d]"%(f,self.firstframetracked,self.lastframetracked,self.firstframebuffered,self.lastframebuffered,self.firstframewritten,self.lastframewritten)
        if f < self.firstframebuffered or f > self.lastframebuffered:

            if DEBUG: print "get_frame: Need to read frame %d from file. buffered frames = [%d,%d], tracked frames = [%d,%d]"%(f,self.firstframebuffered,self.lastframebuffered,self.firstframetracked,self.lastframetracked)

            # read into buffer
            i = int(num.floor( (f - self.firstframetracked) / self.lookupinterval ))

            if DEBUG: print "i = " + str(i)
            if DEBUG: print "lookup element %d = %x"%(i,self.lookup[i])
            self.file.seek(self.lookup[i])
            g = self.firstframetracked + self.lookupinterval*i
            self.buffer = []
            self.firstframebuffered = max(g,f-self.nbuffer+1)
            self.lastframebuffered = min(max(f,g+self.nbuffer-1),self.lastframetracked)
            self.n = self.lastframebuffered-self.firstframebuffered+1
            for h in range(g,self.firstframebuffered):
                s = self.file.readline()
            # add to buffer if necessary
            for h in range(self.firstframebuffered,self.lastframebuffered+1):
                self.buffer.append(self.read_ellipses())

        return self.buffer[f-self.firstframebuffered]

    def get_frames(self,f1,f2):

        # refill buffer if necessary
        ells = []
        for f in range(f1,f2+1):
            ells.append(self.get_frame(f))

        return ells

    def replace_frame(self,ellipses,t):

        if DEBUG: print "Replacing frame %d, frames tracked = [%d,%d], frames buffered = [%d,%d], frames written = [%d,%d]"%(self.firstframetracked,self.lastframetracked,self.firstframebuffered,self.lastframebuffered,self.firstframewritten,self.lastframewritten)

        bufferoff = t - self.firstframebuffered
        if DEBUG: print "bufferoff = " + str(bufferoff)

        # sanity checks
        if t <= self.lastframewritten:
            raise NotImplementedError('Cannot change a frame that is already written to disk')
        if bufferoff < 0 or bufferoff > self.n:
            raise NotImplementedError('Buffer must contain frame to be changed')
        self.buffer[bufferoff] = ellipses

    def append(self,ellipses):
        self.add_frame(ellipses)

    def add_frame(self,ellipses):

        if DEBUG: print "add_frame: %d, initially framestracked = [%d,%d], framesbuffered = [%d,%d], frameswritten = [%d,%d]"%(self.lastframetracked+1,self.firstframetracked,self.lastframetracked,self.firstframebuffered,self.lastframebuffered,self.firstframewritten,self.lastframewritten)

        # write to file:
        # we will be adding frame n+1, so we can write
        # frame n - maxlookback
        if self.n >= self.maxlookback:
            self.lastframewritten+=1
            if DEBUG: print "Writing to file. lastframewritten is now %d"%self.lastframewritten
            if self.islookupframe(self.lastframewritten):
                self.lookup.append(self.file.tell())
                if DEBUG: print "Adding to lookup frame %d (%x)"%(self.lastframewritten,self.lookup[-1])
            if DEBUG: print "Writing frame %d at %x: "%(self.lastframebuffered-self.maxlookback,self.file.tell()) + str(self.buffer[self.n-self.maxlookback])
            self.write_ellipses(self.buffer[self.n-self.maxlookback])

        # update buffer
        if self.n >= self.nbuffer: # should never be >
            tmp = self.buffer.pop(0)
            if DEBUG: print "popping frame %d from buffer: "%self.firstframebuffered + str(tmp)
            self.firstframebuffered+=1
        else:
            self.n += 1

        self.buffer.append(ellipses)
        if DEBUG: print "after adding buffer length = %d, n = %d"%(len(self.buffer),self.n)
        self.lastframebuffered += 1
        self.lastframetracked += 1
        self.nframestracked += 1
        if DEBUG: print "Added to end of buffer[%d] (frame = %d), nflies = %d: "%(len(self.buffer)-1,self.lastframebuffered,len(ellipses)) + str(ellipses)
        if DEBUG: print "framestracked: [%d,%d], frameswritten: [%d,%d], framesbuffered: [%d,%d], nframestracked: %d, n: %d"%(self.firstframetracked,self.lastframetracked,self.firstframewritten,self.lastframewritten,self.firstframebuffered,self.lastframebuffered,self.nframestracked,self.n)

        self.isdatawritten = True

    def __getitem__( self, i ):
        #if DEBUG: print "Calling __getitem__(%d)"%i
        if i < 0:
            j = self.lastframetracked + 1 + i
        else:
            j = self.firstframetracked + i
        return self.get_frame(j)

    def __setitem__( self, i, val ):
        #if DEBUG: print "Calling __setitem__(%d) = "%i + str(val)
        if i < 0:
            j = self.lastframetracked + 1 + i
        else:
            j = self.firstframetracked + i
        if j == self.lastframetracked + 1:
            self.add_frame(val)
        elif j > self.lastframetracked:
            raise NotImplementedError("Trying to set frame %d, lastframetracked = %d")%(j,self.lastframetracked)
        else:
            self.replace_frame(val, j)

    def __len__( self ):
        #if DEBUG: print "calling __len__"
        return self.nframestracked

    def finish_writing(self):
        
        if DEBUG: print "Finishing writing"
        # seek to the end of the file
        self.file.flush()
        # this seek seemed to break things
        #self.file.seek(0,2)

        for f in range(self.lastframewritten+1,self.lastframetracked+1):
            if self.islookupframe(f):
                self.file.flush()
                self.lookup.append(self.file.tell())
            if DEBUG: print "Writing frame %d at %x: "%(f,self.file.tell()) + str(self.buffer[f-self.firstframebuffered])
            self.write_ellipses(self.buffer[f-self.firstframebuffered])
        self.lastframewritten = self.lastframetracked

        self.file.flush()
        self.file.truncate()

        self.isdatawritten = True

    def lookupfloor(self,f,lookupinterval=None,firstframetracked=None):
        if lookupinterval is None:
            lookupinterval = self.lookupinterval
        if firstframetracked is None:
            firstframetracked = self.firstframetracked

        i = int(num.floor( (f - firstframetracked) / lookupinterval ))
        g = firstframetracked + lookupinterval*i
        return (i,g)

    def IsAnnData(self):
        return self.lastframetracked >= self.firstframetracked

    def __del__( self ):
        if hasattr(self,'file') and (self.file is not None) and \
                (not self.file.closed): 
            self.file.close()

    def write_ellipses( self, ellipse_list, fileout=None ):
        """Write one frame of data to already-open file."""

        if fileout is None:
            fileout = self.file

        for ellipse in ellipse_list.itervalues():
            fileout.write( '%f\t%f\t%f\t%f\t%f\t%d\t'%(ellipse.center.x,
                                                       ellipse.center.y,
                                                       ellipse.size.width,
                                                       ellipse.size.height,
                                                       ellipse.angle,
                                                       ellipse.identity) )
        fileout.write( "\n" )

    def write_ellipses_string( self, ellipse_list):
        """Write one frame of data to string."""

        s = ''

        for ellipse in ellipse_list.itervalues():
            s += '%f\t%f\t%f\t%f\t%f\t%d\t'%(ellipse.center.x,
                                             ellipse.center.y,
                                             ellipse.size.width,
                                             ellipse.size.height,
                                             ellipse.angle,
                                             ellipse.identity)
        s += "\n"
        return s

    def WriteAnnHeader( self, start_frame ):
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
        if hasattr(self.bg_imgs,'med'):
            sz = self.bg_imgs.med.size*SIZEOFDOUBLE
            self.file.write("background median:%d\n"%sz)
            self.file.write(self.bg_imgs.med)
        if hasattr(self.bg_imgs,'mean'):
            sz = self.bg_imgs.mean.size*SIZEOFDOUBLE
            self.file.write("background mean:%d\n"%sz)
            self.file.write(self.bg_imgs.mean)
        self.file.write('bg_norm_type:%d\n'%params.bg_norm_type)
        if hasattr(self.bg_imgs,'mad'):
            sz = self.bg_imgs.mad.size*SIZEOFDOUBLE
            self.file.write("background mad:%d\n"%sz)
            self.file.write(self.bg_imgs.mad)
        if hasattr(self.bg_imgs,'std'):
            sz = self.bg_imgs.std.size*SIZEOFDOUBLE
            self.file.write("background std:%d\n"%sz)
            self.file.write(self.bg_imgs.std)
        if hasattr(self.bg_imgs,'hfnorm'):
            sz = self.bg_imgs.hfnorm.size*SIZEOFDOUBLE
            self.file.write("hfnorm:%d\n"%sz)
            self.file.write(self.bg_imgs.hfnorm)
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
        
        self.file.write('data format:%s\n'%dataformatstring)

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

    def ReadAnnHeader( self, doreadbgmodel=True,
                       doreadstartframe=True ):
        """Read header from an annotation file."""

        self.file.seek(0,0)

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
                if DEBUG: print "Found end of annotation header at %x"%self.start_data
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
                params.do_set_circular_arena = bool(int(value))
            elif parameter == 'do_use_morphology':
                params.do_use_morphology = bool(int(value))
            elif parameter == 'opening_radius':
                params.opening_radius = int(value)
                self.bg_imgs.opening_struct = self.bg_imgs.create_morph_struct(params.opening_radius)
            elif parameter == 'closing_radius':
                params.closing_radius = int(value)
                self.bg_imgs.closing_struct = self.bg_imgs.create_morph_struct(params.closing_radius)
            elif parameter == 'bg_algorithm':
                if value == 'median':
                    params.use_median = True
                else:
                    params.use_median = False
            elif parameter == 'background median':
                sz = int(value)
                if doreadbgmodel:
                    self.bg_imgs.med = num.fromstring(self.file.read(sz),'<d')
                    self.bg_imgs.med.shape = params.movie_size
                else:
                    self.file.seek(sz,1)
            elif parameter == 'background mean':
                sz = int(value)
                if doreadbgmodel:
                    self.bg_imgs.mean = num.fromstring(self.file.read(sz),'<d')
                    self.bg_imgs.mean.shape = params.movie_size
                else:
                    self.file.seek(sz,1)
            elif parameter == 'bg_norm_type':
                params.bg_norm_type = int(value)
            elif parameter == 'background mad':
                sz = int(value)
                if doreadbgmodel:
                    self.bg_imgs.mad = num.fromstring(self.file.read(sz),'<d')
                    self.bg_imgs.mad.shape = params.movie_size
                else:
                    self.file.seek(sz,1)
            elif parameter == 'background std':
                sz = int(value)
                if doreadbgmodel:
                    self.bg_imgs.std = num.fromstring(self.file.read(sz),'<d')
                    self.bg_imgs.std.shape = params.movie_size
                else:
                    self.file.seek(sz,1)
            elif parameter == 'hfnorm':
                sz = int(value)
                if doreadbgmodel:
                    self.bg_imgs.hfnorm = num.fromstring(self.file.read(sz),'<d')
                    self.bg_imgs.hfnorm.shape = params.movie_size
                else:
                    self.file.seek(sz,1)
            elif parameter == 'roipolygons':
                sz = int(value)
                if doreadbgmodel:
                    params.roipolygons = pickle.loads(self.file.read(sz))
                else:
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
            elif parameter == 'do_fix_split':
                params.do_fix_split = bool(int(value))
            elif parameter == 'splitdetection_length':
                params.splitdetection_length = int(value)
            elif parameter == 'splitdetection_cost':
                params.splitdetection_cost = float(value)
            elif parameter == 'do_fix_merged':
                params.do_fix_merged = bool(int(value))
            elif parameter == 'mergeddetection_length':
                params.mergeddetection_length = int(value)
            elif parameter == 'mergeddetection_distance':
                params.mergeddetection_distance = float(value)
            elif parameter == 'do_fix_spurious':
                params.do_fix_spurious = bool(int(value))
            elif parameter == 'spuriousdetection_length':
                params.spuriousdetection_length = int(value)
            elif parameter == 'do_fix_lost':
                params.do_fix_lost = bool(int(value))
            elif parameter == 'lostdetection_length':
                params.lostdetection_length = int(value)
            elif parameter == 'movie_name':
                params.annotation_movie_name = value
            elif parameter == 'start_frame':
                if doreadstartframe:
                    params.start_frame = int(value)
            elif parameter == 'data format':
                self.n_fields = len(value.split())

        if doreadbgmodel:
            self.bg_imgs.updateParameters()

    def ReadSettings( self ):
        """Read header from an annotation file."""

        self.ReadAnnHeader(doreadbgmodel=False,doreadstartframe=False)

    def ReplaceIds(self,ells):
        new_ells = ell.TargetList()

        # count, replace identities
        for (i,e) in ells.iteritems():
            if i not in self.idtable:
                if DEBUG: print "identity %d is not in idtable: "%i + str(self.idtable)
                self.idtable[i] = params.nids
                params.nids += 1
                if DEBUG: print "incrementing params.nids to %d"%params.nids
            e.identity = self.idtable[i]
            new_ells.append(e)

        return new_ells

    def CountIds(self,ells):

        # count, replace identities
        for (i,e) in ells.iteritems():
            if i not in self.idtable:
                if DEBUG: print "identity %d is not in idtable: "%i + str(self.idtable)
                self.idtable[i] = params.nids
                params.nids += 1
                if DEBUG: print "incrementing params.nids to %d"%params.nids

    def ParseData( self, line, doreplaceids=False):
        """Split a line of annotation data into per-fly information.
        Returns an EllipseList. Allows passing in an EllipseList, which is
        overwritten and returned, to avoid memory reallocation."""
        
        fly_sep = line.split()
        if (len(fly_sep) % self.n_fields) > 0:
            print "Error reading trajectories from annotation file. "
            print "line = %s"%line
            print "parsed as " + str(fly_sep)
            raise Exception, 'Error reading trajectories from annotation file'

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

            ellipses.append( new_ellipse )

        return ellipses
    
    def read_ellipses_string(self,s):
        return self.ParseData(s)
    
    def read_ellipses(self,filein=None):
        if filein is None:
            filein = self.file
        return self.ParseData(filein.readline())

    # need to fix this
    def WriteMAT( self, filename, dosavestamps=False ):
        """Writes a MATLAB .mat file from the data sent in. Expects a
        list of EllipseLists."""
        # how many targets per frame
        nframes = self.lastframetracked - self.firstframetracked + 1
        startframe = self.firstframetracked
        ntargets = num.ones( nframes )

        # we are going to store everything in memory. this may be prohibitive
        if self.nbuffer < nframes:
            msgtxt = 'Saving to a mat file will require loading all %d frames into memory, which may require a prohibitive amount of RAM. '%nframes
            if params.interactive and not params.batch_executing:
                a = wx.MessageBox(msgtxt + 'Proceed?',"Store to MAT file?",
                                  wx.YES|wx.NO)
                if a == wx.NO:
                    return
            else:
                print msgtxt

        # first get ntargets to allocate other matrices
        off = self.firstframetracked
        for i in range(nframes):
            ntargets[i] = len(self.get_frame(off+i))
        z = num.sum(ntargets)

        # allocate arrays as all NaNs
        x_pos = num.ones( z ) * num.nan
        y_pos = num.ones( z ) * num.nan
        maj_ax = num.ones( z ) * num.nan
        min_ax = num.ones( z ) * num.nan
        angle = num.ones( z ) * num.nan
        identity = num.ones( z ) * num.nan

        # store data in arrays
        # go through all the frames
        i = 0
        for j in range( nframes ):
            ells = self.get_frame(off+j)
            for ee in ells.itervalues():
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
            if DEBUG: print "reading in movie to save stamps"
            stamps = params.movie.get_some_timestamps(t1=startframe,t2=startframe+nframes)
            save_dict['timestamps'] = stamps

        scipy.io.savemat( filename, save_dict, oned_as='column' )
        # menu_file_write_timestamps
    
    def tell():
        return self.file.tell()

    def get_frame_prev(self):
        return self.get_frame(self.lastframetracked-1)
    def get_frame_prevprev(self):
        return self.get_frame(self.lastframetracked-2)

    def GetNewId(self):
        if len(self.recycledids) > 0:
            newid = self.recycledids.pop()
            if DEBUG: print "Recycled id %d"%newid
        else:
            newid = params.nids
            if DEBUG: print "Used new id %d"%newid
            params.nids+=1
        return newid

    def RecycleId(self,id):
        self.recycledids.append(id)
        if DEBUG: print "Recycling id %d"%id

    def rename_file(self,newfilename=None):
        oldfile = self.file
        oldfilename = self.filename
        if newfilename is None:
            self.file = tempfile.NamedTemporaryFile(suffix='.ann')
            self.filename = self.file.name
        else:
            self.filename = newfilename
            self.file = open(self.filename,'wb+')
        if DEBUG: print "rename_file: created new file %s"%self.filename
        # seek to the start
        oldfile.seek(0,0)
        # copy
        shutil.copyfileobj(oldfile,self.file)
        oldfile.close()

def LoadSettings(filename, bg_imgs, doreadbgmodel=False):
    tmpannfile = AnnotationFile(filename,
                                bg_imgs,
                                justreadheader=True,
                                doreadbgmodel=doreadbgmodel)
