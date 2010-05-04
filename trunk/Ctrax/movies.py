# movies.py
# KMB 11/06/2008

import chunk
import numpy as num
import struct
#import threading
import wx
import os
from params import params
from draw import annotate_bmp
import pyglet.media as media
from version import DEBUG_MOVIES as DEBUG

# version of sbfmf for writing
__version__ = "0.3b"

#import static_background_fmf as sbfmf
#import motmot.FlyMovieFormat.FlyMovieFormat as fmf
import FlyMovieFormat as fmf
try:
    from FlyMovieFormat import NoMoreFramesException
except ImportError:
    class NoMoreFramesException (Exception): pass
try:
    import motmot.ufmf.ufmf as ufmf
except ImportError:
    pass

class Movie:
    """Generic interface for all supported movie types.
    Also makes reading movies thread-safe."""
    def __init__( self, filename, interactive ):
        """Figure out file type and initialize reader."""
        self.interactive = interactive

        #self.lock = threading.Lock()

        # read FlyMovieFormat
        (tmp,ext) = os.path.splitext(filename)
        if ext == '.fmf':
            self.type = 'fmf'
            try:
                self.h_mov = fmf.FlyMovie( filename )
            except NameError:
                if self.interactive:
                    wx.MessageBox( "Couldn't open \"%s\"\n(maybe FMF is not installed?)"%(filename), "Error", wx.ICON_ERROR )
                raise
            except IOError:
                if self.interactive:
                    wx.MessageBox( "I/O error opening \"%s\""%(filename), "Error", wx.ICON_ERROR )
                raise
        elif ext == '.sbfmf':
            self.type = 'sbfmf'
            try:
                self.h_mov = fmf.FlyMovie(filename)
            except NameError:
                if self.interactive:
                    wx.MessageBox( "Couldn't open \"%s\"\n(maybe FMF is not installed?)"%(filename), "Error", wx.ICON_ERROR )
                raise
            except IOError:
                if self.interactive:
                    wx.MessageBox( "I/O error opening \"%s\""%(filename), "Error", wx.ICON_ERROR )
                raise
        elif ext == '.ufmf':
            self.type = 'ufmf'
            try:
                self.h_mov = ufmf.FlyMovieEmulator(filename)
            except NameError:
                if self.interactive:
                    wx.MessageBox( "Couldn't open \"%s\"\n(maybe UFMF is not installed?)"%(filename), "Error", wx.ICON_ERROR )
                raise
            except IOError:
                if self.interactive:
                    wx.MessageBox( "I/O error opening \"%s\""%(filename), "Error", wx.ICON_ERROR )
                raise
            except ufmf.ShortUfmfFileError:
                if self.interactive:
                    wx.MessageBox( "Error opening \"%s\". Short ufmf file."%(filename), "Error", wx.ICON_ERROR )
                raise
            except ufmf.CorruptIndexError:
                if self.interactive:
                    wx.MessageBox( "Error opening \"%s\". Corrupt file index."%(filename), "Error", wx.ICON_ERROR )
                raise
            except ufmf.InvalidMovieFileException:
                if self.interactive:
                    wx.MessageBox( "Error opening \"%s\". Invalid movie file."%(filename), "Error", wx.ICON_ERROR )
                raise
        # read AVI
        elif ext == '.avi':
            try:
                self.h_mov = Avi( filename, fmfmode=True )
                self.type = 'avi'
            except:
                try:
                    self.h_mov = CompressedAvi( filename )
                    self.type = 'avbin'
                except:
                    if self.interactive:
                        wx.MessageBox( "Failed opening file \"%s\"."%(filename), "Error", wx.ICON_ERROR )
                    raise
                else:
                    if params.interactive:
                        wx.MessageBox("Your movie is most likely compressed. Out-of-order frame access (e.g. dragging the frame slider toolbars around) will be slow. At this time, the frame chosen to be displayed may be off by one or two frames, i.e. may not line up perfectly with computed trajectories.","Warning",wx.ICON_WARNING)
            if params.interactive and self.h_mov.bits_per_pixel == 24:
                wx.MessageBox( "Currently, RGB movies are immediately converted to grayscale. All color information is ignored.", "Warning", wx.ICON_WARNING )

        # unknown movie type
        else:
            try:
                self.h_mov = CompressedAvi( filename )
                self.type = 'avbin'
            except:
                if self.interactive:
                    wx.MessageBox( "Failed opening file \"%s\"."%(filename), "Error", wx.ICON_ERROR )
                raise
            else:
                if params.interactive:
                    wx.MessageBox("Your movie is most likely compressed. Out-of-order frame access (e.g. dragging the frame slider toolbars around) will be slow. At this time, the frame chosen to be displayed may be off by one or two frames, i.e. may not line up perfectly with computed trajectories.","Warning",wx.ICON_WARNING)
        # add a buffer of the current frame
        self.bufferedframe_im = None
        self.bufferedframe_stamp = None
        self.bufferedframe_num = None

    def get_frame( self, framenumber ):
        """Return numpy array containing frame data."""
        # check to see if we have buffered this frame
        if framenumber == self.bufferedframe_num:
            return (self.bufferedframe_im.copy(),self.bufferedframe_stamp)

        #self.lock.acquire()
        try:
            try: frame, stamp = self.h_mov.get_frame( framenumber )
            except (IndexError, NoMoreFramesException):
                if self.interactive:
                    wx.MessageBox( "Frame number %d out of range"%(framenumber), "Error", wx.ICON_ERROR )
                raise
            except (ValueError, AssertionError):
                if self.interactive:
                    wx.MessageBox( "Error reading frame %d"%(framenumber), "Error", wx.ICON_ERROR )
                raise
            else:
                # store the current frame in the buffer
                self.bufferedframe_im = frame.copy()
                self.bufferedframe_stamp = stamp
                self.bufferedframe_num = framenumber

                return frame, stamp
        finally:
            pass
            #self.lock.release()

    def get_n_frames( self ): return self.h_mov.get_n_frames()
    def get_width( self ): return self.h_mov.get_width()
    def get_height( self ): return self.h_mov.get_height()
    def get_some_timestamps( self, t1=0, t2=num.inf ):
        t2 = min(t2,self.get_n_frames())
        timestamps = self.h_mov.get_all_timestamps()
        timestamps = timestamps[t1:t2]
        return timestamps

    def writesbfmf_start(self,bg,filename):

        # write from start_frame to nframes-1
        self.nframescompress = self.get_n_frames() - params.start_frame

        # allocate array to store the addresses of each frame
        self.writesbfmf_framestarts = num.zeros(self.nframescompress)

        # open the output file
        self.writesbfmf_outfilename = filename
        self.outfile = open(self.writesbfmf_outfilename,'wb')

        # write the header
        self.writesbfmf_writeheader(bg)

    def writesbfmf_isopen(self):
        if self.outfile is None:
            return False
        return (not self.outfile.closed)

    def writesbfmf_restart(self,frame,bg,filename):

        self.outfile = None

        self.writesbfmf_outfilename = filename

        self.nframescompress = self.get_n_frames() - params.start_frame

        # allocate array to store the addresses of each frame
        self.writesbfmf_framestarts = num.zeros(self.nframescompress)

        # move the file to a temporary file
        tmpfilename = 'tmp_ctrax_writesbfmf.sbfmf'
        os.rename(filename,tmpfilename)

        # open the old file for reading
        inmovie = Movie( tmpfilename, self.interactive )
            
        # open the output file
        self.outfile = open(filename,"wb")
        
        # rewrite the header
        self.writesbfmf_writeheader(bg)

        # last frame to copy
        i = frame - params.start_frame - 1
        firstaddr = inmovie.h_mov.framelocs[0]
        lastaddr = inmovie.h_mov.framelocs[i]

        self.writesbfmf_framestarts[:i+1] = inmovie.h_mov.framelocs[:i+1]

        if DEBUG: print "copied framestarts: "
        if DEBUG: print str(self.writesbfmf_framestarts[:i+1])

        # seek to the first frame
        inmovie.h_mov.seek(0)

        # copy in pages of size pagesize
        pagesize = int(2**20)
        for j in range(firstaddr,lastaddr,pagesize):
            if DEBUG: print "writing page at %d"%inmovie.h_mov.file.tell()
            buf = inmovie.h_mov.read_some_bytes(pagesize)
            self.outfile.write(buf)

        # write last page
        if j < lastaddr:
            if DEBUG: print "writing page at %d"%inmovie.h_mov.file.tell()
            buf = inmovie.h_mov.read_some_bytes(int(lastaddr-pagesize))
            self.outfile.write(buf)

        # close the input movie and delete
        inmovie.h_mov.close()
        os.remove(tmpfilename)

    def writesbfmf_close(self,frame):

        # write the index
        self.writesbfmf_writeindex(frame)

        # close the file
        self.outfile.close()

    def writesbfmf_writeindex(self,frame):
        """
        Writes the index at the end of the file. Index consists of nframes unsigned long longs (Q),
        indicating the positions of each frame
        """
        # write the index
        indexloc = self.outfile.tell()
        nframeswrite = frame - params.start_frame + 1
        if DEBUG: print "writing index, nframeswrite = %d"%nframeswrite
        for i in range(nframeswrite):
            self.outfile.write(struct.pack("<Q",self.writesbfmf_framestarts[i]))

        # write the location of the index
        self.outfile.seek(self.writesbfmf_indexptrloc)
        self.outfile.write(struct.pack("<Q",indexloc))            

        # write the number of frames
        self.outfile.seek(self.writesbfmf_nframesloc)
        self.outfile.write(struct.pack("<I",nframeswrite))


    def writesbfmf_writeheader(self,bg):
        """
        Writes the header for the file. Format:
        Number of bytes in version string: (I = unsigned int)
        Version Number (string of specified length)
        Number of rows (I = unsigned int)
        Number of columns (I = unsigned int)
        Number of frames (I = unsigned int)
        Difference mode (I = unsigned int):
          0 if light flies on dark background, unsigned mode
          1 if dark flies on light background, unsigned mode
          2 if other, signed mode
        Location of index (Q = unsigned long long)
        Background image (ncols * nrows * double)
        Standard deviation image (ncols * nrows * double)
        """

        self.nr = self.get_height()
        self.nc = self.get_width()

        # write the number of columns, rows, frames, difference mode
        if params.bg_type == params.BG_TYPE_LIGHTONDARK:
            difference_mode = 0
        elif params.bg_type == params.BG_TYPE_DARKONLIGHT:
            difference_mode = 1
        else:
            difference_mode = 2
        self.outfile.write(struct.pack("<I",len(__version__)))
        self.outfile.write(__version__)
        self.outfile.write(struct.pack("<2I",int(self.nr),int(self.nc)))
        self.writesbfmf_nframesloc = self.outfile.tell()
        self.outfile.write(struct.pack("<2I",int(self.nframescompress),int(difference_mode)))

        if DEBUG: print "writeheader: nframescompress = " + str(self.nframescompress)

        # compute the location of the standard deviation image
        stdloc = self.outfile.tell() + struct.calcsize("B")*self.nr*self.nc

        # compute the location of the first frame
        ffloc = stdloc + struct.calcsize("d")*self.nr*self.nc

        # where do we write the location of the index -- this is always the same
        self.writesbfmf_indexptrloc = self.outfile.tell()

        # write a placeholder for the index location
        self.outfile.write(struct.pack("<Q",0))

        # write the background image
        self.outfile.write(bg.center)

        # write the standard deviation image
        self.outfile.write(bg.dev)

    def writesbfmf_writeframe(self,isfore,im,stamp,currframe):

        if DEBUG: print "writing frame %d"%currframe

        tmp = isfore.copy()
        tmp.shape = (self.nr*self.nc,)
        i, = num.nonzero(tmp)

        # values at foreground pixels
        v = im[isfore]

        # number of foreground pixels
        n = len(i)

        # store the start of this frame
        j = currframe - params.start_frame
        self.writesbfmf_framestarts[j] = self.outfile.tell()

        if DEBUG: print "stored in framestarts[%d]"%j

        # write number of pixels and time stamp
        self.outfile.write(struct.pack("<Id",n,stamp))

        i = i.astype(num.uint32)

        self.outfile.write(i)
        self.outfile.write(v)


"""
AVI class; written by JB and KMB, altered by Don Olbris.

Don's changes:
important alterations from version I received:
- allows fccHandler = "" or all nulls
- allows width to vary to match actual frame size (suspect that avi pads to 
    multiples of four?) (and let get_width return actual width)
- added various derived attributes to better masquerade as an FMF file
- added an "fmf" mode; in this mode, it reshapes array to same shape as fmf
    (row/column order swapped)
- added a seek() method

- I want to make this width change more transparent, but I keep
    running into related shape issues

- Avi class is still byte-order undefined, but at least it's read-only

"""

class Avi:
    """Read uncompressed AVI movies."""
    def __init__( self, filename, fmfmode=False ):
        
        self.fmfmode = fmfmode
        self.issbfmf = False
 
        # need to open in binary mode to support Windows:
        self.file = open( filename, 'rb' )
        try:
            self.read_header()
        except Exception, details:
            if DEBUG: print( "error reading uncompressed AVI:" )
            if DEBUG: print( details )
            raise
        
        # added to help masquerade as FMF file:
        self.filename = filename
        self.chunk_start = self.data_start
        # this is a mystery, but I think it's 8 for avi: seems to be offset
        #   from beginning of "chunk" to beginning of array data within chunk
        self.timestamp_len = 8
        if hasattr(self, "newwidth"):
            self.bytes_per_chunk = (self.height * self.newwidth) + self.timestamp_len
        else:
            self.bytes_per_chunk = self.buf_size + self.timestamp_len
        #self.bits_per_pixel = 8 

        self._all_timestamps = None
    
    def get_all_timestamps(self):
        if self._all_timestamps is None:
            self._all_timestamps = []
            for f in range(self.n_frames):
                try:
                    (frame,ts) = self.get_frame(f)
                except:
                    continue
                self._all_timestamps.append(ts)
        return self._all_timestamps

    def read_header( self ):

        # read RIFF then riffsize
        RIFF, riff_size, AVI = struct.unpack( '4sI4s', self.file.read( 12 ) )
        if not RIFF == 'RIFF':
            raise TypeError("Invalid AVI file. Must be a RIFF file.")
        if not AVI == 'AVI ':
            raise TypeError("Invalid AVI file. File type must be AVI .")

        # read hdrl
        LIST, hdrl_size, hdrl = \
            struct.unpack( '4sI4s', self.file.read( 12 ) )
        hdrlstart = self.file.tell() - 4

        if not LIST == 'LIST' or not hdrl == 'hdrl':
            raise TypeError("Invalid AVI file. Did not find header list. %s should = LIST and %s should = hdrl"%(LIST,hdrl))

        # read avih 
        avih, avih_size = struct.unpack( '4sI', self.file.read( 8 ) )
        if not avih == 'avih':
            raise TypeError("Invalid AVI file. Did not find avi header.")
        avihchunkstart = self.file.tell()

        # read microsecperframe
        self.frame_delay_us, = struct.unpack('I',self.file.read(4))
        
        # skip to nframes
        self.file.seek(3*4,1)
        self.n_frames, = struct.unpack('I',self.file.read(4))

        # skip to width
        self.file.seek(3*4,1)
        self.width,self.height = struct.unpack('2I',self.file.read(8))

        if DEBUG: print "width = %d, height = %d"%(self.width,self.height)
        if DEBUG: print "n_frames = %d"%self.n_frames

        # skip the rest of the aviheader
        self.file.seek(avihchunkstart+avih_size,0)

        LIST, stream_listsize, strl = \
            struct.unpack( '4sI4s', self.file.read( 12 ) )

        if not LIST == 'LIST' or not strl == 'strl':
            raise TypeError("Invalid AVI file. Did not find stream list.")

        strh, strh_size = struct.unpack( '4sI', self.file.read( 8 ) )
        if not strh == 'strh':
            raise TypeError("Invalid AVI file. Did not find stream header.")

        strhstart = self.file.tell()

        # read stream type, fcc handler
        vids, fcc = struct.unpack( '4s4s', self.file.read( 8 ) )
        # check for vidstream
        if not vids == 'vids':
            raise TypeError("Unsupported AVI file type. First stream found is not a video stream.")
        # check fcc
        if fcc not in ['DIB ', '\x00\x00\x00\x00', "", "RAW ", "NONE", chr(24)+"BGR"]:
            raise TypeError("Unsupported AVI file type %s, only uncompressed AVIs supported."%fcc)

        # skip the rest of the stream header
        self.file.seek(strhstart+strh_size,0)
        
        strf, strf_size = struct.unpack( '4sI', self.file.read( 8 ) )
        if not strf == "strf":
            raise TypeError("Invalid AVI file. Did not find strf.")

        strfstart = self.file.tell()
        bitmapheadersize, = struct.unpack('I',self.file.read(4))

        # skip width, height, planes
        self.file.seek(4*2+2,1)

        # read in bits per pixel
        self.bits_per_pixel, = struct.unpack('H',self.file.read(2))
        if DEBUG: print "bits_per_pixel = %d"%self.bits_per_pixel

        # is this an indexed avi?
        colormapsize = (strf_size - bitmapheadersize)/4
        if colormapsize > 0:
            self.isindexed = True
            self.file.seek(strfstart+bitmapheadersize,0)
            self.colormap = num.frombuffer(self.file.read(4*colormapsize),num.uint8)
            self.colormap = self.colormap.reshape((colormapsize,4))
            self.colormap = self.colormap[:,:-1]
        else:
            self.isindexed = False

        # skip the rest of the strf
        self.file.seek(hdrlstart+hdrl_size,0)

        #if self.isindexed:
        #    raise TypeError("Unsupported AVI file type: indexed colormap. Only RGB and grayscale, uncompressed AVIs are currently supported")
        
        while True:
            # find LIST chunk
            LIST,movilist_size = struct.unpack( '4sI', self.file.read( 8 ) )
            if LIST == 'LIST':
                # find movi
                movistart = self.file.tell()
                movi, = struct.unpack('4s',self.file.read(4))
                if movi == 'movi':
                    break
                else:
                    self.file.seek(-4,1)
            # found some other chunk, seek past
            self.file.seek(movilist_size,1)            

        # sometimes there is junk here
        #if LIST == 'JUNK':
        #    self.file.seek(movilist_size,1)
        #    LIST,movilist_size = struct.unpack( '4sI', self.file.read( 8 ) )
        #movistart = self.file.tell()

        #if not LIST == 'LIST':
        #    raise TypeError("Invalid AVI file. Did not find movie LIST.")

        #movi, = struct.unpack('4s',self.file.read(4))

        if not movi == 'movi':
            raise TypeError("Invalid AVI file. Did not find movi, found %s."%movi)

        self.data_start = self.file.tell()

        # read one frame's header to check
        this_frame_id, frame_size = struct.unpack( '4sI', self.file.read( 8 ) )
        if DEBUG: print "frame_size = " + str(frame_size)
        depth = self.bits_per_pixel/8

        # figure out padding
        unpaddedframesize = self.width*self.height*depth
        if unpaddedframesize == frame_size:
            # no padding
            self.padwidth = 0
            self.padheight = 0
        elif unpaddedframesize + self.width*depth == frame_size:
            self.padwidth = 0
            self.padheight = 1
        elif unpaddedframesize + self.height*depth == frame_size:
            self.padwidth = 1
            self.padheight = 0
        else:
            raise TypeError("Invalid AVI file. Frame size (%d) does not match width * height * bytesperpixel (%d*%d*%d)."%(frame_size,self.width,self.height,depth))

        if self.isindexed:
            self.format = 'INDEXED'
        elif self.bits_per_pixel == 8:
            self.format = 'MONO8'
        elif self.bits_per_pixel == 24:
            self.format = 'RGB'
        else:
            raise TypeError("Unsupported AVI type. bitsperpixel must be 8 or 24, not %d."%self.bits_per_pixel)

        if DEBUG: print "format = " + str(self.format)

        # set buf size
        self.buf_size = frame_size
        #self.buf_size = self.width*self.height*self.bits_per_pixel/8


    def old_read_header( self ):

        if DEBUG: print "reading chunk header"
        # RIFF header
        file_type, riff_size = struct.unpack( '4sI', self.file.read( 8 ) )
        assert file_type == 'RIFF'
        stream_type = self.file.read( 4 )
        assert stream_type == 'AVI '
        header_list, header_listsize, header_listtype = \
                     struct.unpack( '4sI4s', self.file.read( 12 ) )
        assert header_list == 'LIST' and header_listtype == 'hdrl'
        #size 4588 (fmf 1222)
        avi_str, avi_note = struct.unpack( '4sI', self.file.read( 8 ) )
        assert avi_str == 'avih'

        if DEBUG: print "1"

        # AVI header
        avi_header = self.file.read( 56 )
        self.frame_delay_us, \
                          AVI_data_rate, \
                          padding_size, \
                          AVI_flags, \
                          self.n_frames, \
                          n_preview_streams, \
                          n_data_streams, \
                          avi_buf_size, \
                          self.width, \
                          self.height, \
                          self.time_scale, \
                          self.data_rate, \
                          self.start_time, \
                          self.AVI_chunk_size \
                          = struct.unpack( '14I', avi_header )
        #10000 100000000 0 16 100 0 1 1310720 1280 1024 100 10000 0 99
        if n_data_streams != 1:
            raise TypeError( "file must contain only one data stream" )
        if avi_buf_size != 0: self.buf_size = avi_buf_size

        if DEBUG: print "2"

        # stream header
        stream_list, stream_listsize, stream_listtype = \
                     struct.unpack( '4sI4s', self.file.read( 12 ) )
        assert stream_list == 'LIST' and stream_listtype == 'strl'
        #size 4244 (fmf 1146)
        stream_str, stream_note = struct.unpack( '4sI', self.file.read( 8 ) )
        assert stream_str == 'strh'
        
        stream_header = self.file.read( 56 )
        fccType, \
                 fccHandler, \
                 stream_flags, \
                 priority, \
                 frames_interleave, \
                 stream_scale, \
                 stream_rate, \
                 stream_start, \
                 stream_length, \
                 stream_buf_size, \
                 stream_quality, \
                 stream_sample_size, \
                 x,y,w,h \
                 = struct.unpack( '4s4s10I4H', stream_header )
        #vids DIB  0 0 0 100 10000 0 99 1310720 100 0 0 0 1280 1024
        #vids DIB  0 0 0 1 50 0 301 235200 0 0 0 0 560 420
        
        if DEBUG: print "3"

        if fccType != 'vids':
            raise TypeError( "stream type must be video" )
        # Reiser's current avi's have fccHandler = '\x00\x00\x00\x00'; I
        #   saw a reference that "" is also acceptable (djo)
        if fccHandler not in ['DIB ', '\x00\x00\x00\x00', ""]:
            print "video must be uncompressed; found fccHandler %s" % fccHandler
            raise ValueError( "video must be uncompressed; found fccHandler %s" % fccHandler)
            # raise TypeError( "video must be uncompressed" )
        if DEBUG: print "stream_buf_size = " + str(stream_buf_size)
        if stream_buf_size != 0:
            if hasattr( self, 'buf_size' ):
                if DEBUG: print "buf_size = %d should = stream_buf_size = %d"%(self.buf_size,stream_buf_size)
                assert self.buf_size == stream_buf_size
            else:
                self.buf_size = stream_buf_size

        # bitmap header

        if DEBUG: print "3.5"

        bmp_str, bmp_note = struct.unpack( '4sI', self.file.read( 8 ) )
        if DEBUG: print "bmp_str = " + str(bmp_str)

        assert bmp_str == 'strf'

        if DEBUG: print "4"

        bmp_header = self.file.read( 40 )
        self.bmp_size, \
                       bmp_width, \
                       bmp_height, \
                       bmp_planes, \
                       bmp_bitcount, \
                       crap, \
                       bmp_buf_size, \
                       xpels_per_meter, \
                       ypels_per_meter, \
                       color_used, \
                       color_important \
                       = struct.unpack( 'I2i2H6i', bmp_header )
        #40 1280 1024 1 8 0 1310720 0 0 256 0
        assert bmp_width == self.width and bmp_height == self.height
        if bmp_buf_size != 0:
            if hasattr( self, 'buf_size' ):
                assert self.buf_size == bmp_buf_size
            else:
                self.buf_size == stream_buf_size
        if not hasattr( self, 'buf_size' ):
            # just a guess -- should be OK if 8-bit grayscale
            self.buf_size = self.height * self.width

        if DEBUG: print "5"

        # skip extra header crap
        movie_list = ''
        movie_listtype = ''
        while movie_list != 'LIST':
            s = ''
            EOF_flag = False
            while s.find( 'movi' ) < 0 and not EOF_flag:
                p = self.file.tell()
                s = self.file.read( 128 )
                '''
                # debug: could 'movi' fall across a the boundary of a 128 byte block?
                if p > 5200:
                    print "failed to find 'movi' as needed"
                    sys.exit()
                print "location:", p
                # print "data:", s
                '''
                if s == '': EOF_flag = True
            if EOF_flag: break
            self.file.seek( p )
            self.file.read( s.find( 'movi' ) - 8 )
            movie_list, movie_listsize, movie_listtype = \
                        struct.unpack( '4sI4s', self.file.read( 12 ) )
        assert movie_list == 'LIST' and movie_listtype == 'movi'

        # beginning of data blocks
        self.data_start = self.file.tell()

        if DEBUG: print "6"
        
        # attempt to do frame-size vs width check here:
        this_frame_id, frame_size = struct.unpack( '4sI', self.file.read( 8 ) )
        if frame_size == self.width * self.height:
            # this is fine
            pass
        elif frame_size == 3 * self.width * self.height:
            # probably RGB, no good:
            raise TypeError( "movie must be grayscale" )
        else:
            # frame size doesn't match; for this exercise, pretend the height is 
            #   right and see if width is integral and within 10 of expected;
            #   if so, use that; otherwise, error (djo)
            # raise ValueError( "frame size %d doesn't make sense: movie must be 8-bit grayscale"%(frame.size) )
            if frame_size % self.height == 0:
                self.newwidth = frame_size / self.height
                if abs(self.newwidth - self.width) > 10:
                    raise ValueError("apparent new width = %d; expected width = %d"
                        % (self.height, self.newwidth))
            else:
                raise ValueError("apparent new width is not integral; mod = %d" % (frame_size % self.height))
        
        if DEBUG: print "8"

        
        # end read_header()

    def get_frame( self, framenumber ):
        """Read frame from file and return as NumPy array."""

        if DEBUG: print "uncompressed get_frame(%d)"%framenumber

        if framenumber < 0: raise IndexError
        if framenumber >= self.n_frames: raise NoMoreFramesException
        
        # read frame from file
        self.file.seek( self.data_start + (self.buf_size+8)*framenumber )
       
        # rest of this function has been moved into get_next_frame(), which
        #   pretty much just reads without the seek
        
        return self.get_next_frame()
    
    def get_next_frame(self):
        """returns next frame"""
       
        currentseekloc = self.file.tell()
         
        this_frame_id, frame_size = struct.unpack( '4sI', self.file.read( 8 ) )

        if frame_size != self.buf_size:
            raise ValueError( "Frame size does not equal buffer size; movie must be uncompressed" )
        if not hasattr( self, 'frame_id' ):
            self.frame_id = this_frame_id
        elif this_frame_id != self.frame_id:
            raise ValueError( "error seeking frame start: unknown data header" )
        frame_data = self.file.read( frame_size )

        # frame id: seems to be recorded from first frame, then compared at
        #   each successive frame; just a marker to be checked?
        
        # make frame into numpy array
        frame = num.fromstring( frame_data, num.uint8 )

        width = self.width + self.padwidth
        height = self.height + self.padheight
        if self.isindexed:
            frame = self.colormap[frame,:]
            frame.resize((width,height,3))
            frame = frame[:self.width,:self.height,:]
            tmp = frame.astype(float)
            tmp = tmp[:,:,0]*.3 + tmp[:,:,1]*.59 + tmp[:,:,2]*.11
            tmp = tmp.T
            frame = tmp.astype(num.uint8)
        elif frame.size == width*height:
            frame.resize( (height, width) )
            frame = frame[:self.height,:self.width]
        elif frame.size == width*height*3:
            frame.resize( (height, width*3) )
            tmp = frame.astype(float)
            tmp = tmp[:,2:width*3:3]*.3 + \
                tmp[:,1:width*3:3]*.59 + \
                tmp[:,0:width*3:3]*.11
            tmp = tmp[:self.height,:self.width]
            frame = tmp.astype(num.uint8)
            #frame = imops.to_mono8( 'RGB24', frame )
            #raise TypeError( "movie must be grayscale" )
        else:
            # frame size doesn't match; for this exercise, pretend the height is 
            #   right and see if width is integral and within 10 of expected;
            #   if so, use that; otherwise, error (djo)
            # raise ValueError( "frame size %d doesn't make sense: movie must be 8-bit grayscale"%(frame.size) )
            if frame.size % height == 0:
                self.newwidth = frame.size / height
                if abs(self.newwidth - width) < 10:
                    frame.resize((self.newwidth, height))
                    frame = frame[:self.width,:self.height]
                else:
                    raise ValueError("apparent new width = %d; expected width = %d"
                        % (height, self.newwidth))
            else:
                raise ValueError("apparent new width is not integral; mod = %d" % (frame.size % height))
            
        # make up a timestamp based on the file's stated framerate
        
        # since we don't know the frame number, back it out from
        #   the file location (cleverly grabbed before we started
        #   reading the frame):
        
        framenumber = (currentseekloc - self.chunk_start) / self.bytes_per_chunk
        
        if self.frame_delay_us != 0:
            stamp = framenumber * self.frame_delay_us / 1e6
        elif self.time_scale != 0:
            stamp = framenumber * self.data_rate / float(self.time_scale)
        else:
            stamp = framenumber / 24
            # should raise warning or error here?
            # this might screw up playback, at least
        
        
        #if self.fmfmode:
        #    # row/column order is swapped:
        #    shape = frame.shape
        #    frame.shape = (shape[1], shape[0])
        
        return frame, stamp
        
        # end get_next_frame()
    
    def get_n_frames( self ): 
        return self.n_frames
    
    def get_width( self ): 
        if hasattr(self, "newwidth"):
            return self.newwidth
        else:
            return self.width
    
    def get_height( self ): 
        return self.height
    
    def seek(self,frame_number):
        if frame_number < 0:
            frame_number = self.n_frames + frame_number
        seek_to = self.chunk_start+self.bytes_per_chunk*frame_number
        self.file.seek(seek_to)
    
    # end class Avi


class CompressedAvi:

    """Use pyglet.media to read compressed avi files"""
    def __init__(self,filename):


        if DEBUG: print 'Trying to read compressed AVI'
        self.issbfmf = False
        self.source = media.load(filename)
        self.ZERO = 1./1000000.
        # for some video types, if we don't do a read before 
        # seeking, bad things happen
        im = self.source.get_next_video_frame()
        ts = self.source.get_next_video_timestamp()
        if DEBUG: print 'first frame read results in im = ' + str(im)
        if DEBUG: print 'and ts = ' + str(ts)
        self.source._seek(self.ZERO)
        self.start_time = self.source.get_next_video_timestamp()
        if DEBUG: print 'First seek to ZERO succeeded, start_time = ' + str(self.start_time)

        # estimate properties of the video that would be nice 
        # to be able to read in -- i have no idea how accurate 
        # these estimates are
        self.duration_seconds = self.source.duration
        if DEBUG: print 'duration_seconds = ' + str(self.duration_seconds)
        self._estimate_fps()
        if DEBUG: print 'fps estimated to be' + str(self.fps)
        self.n_frames = int(num.floor(self.duration_seconds * self.fps))
        if DEBUG: print "n_frames estimated to be %d"%self.n_frames
        self.frame_delay_us = 1e6 / self.fps
        self._estimate_keyframe_period()
        if DEBUG: print 'keyframe period estimate to be ' + str(self.keyframe_period_s) + ' s, ' + str(self.keyframe_period) + ' frames'
        
        # added to help masquerade as FMF file:
        self.filename = filename

        # read in the width and height of each frame
        self.width = self.source.video_format.width
        self.height = self.source.video_format.height

        self.MAXBUFFERSIZE = num.round(200*1000*1000/self.width/self.height)
        self.buffersize = min(self.MAXBUFFERSIZE,self.keyframe_period)
        if DEBUG: print 'buffersize set to ' + str(self.buffersize)

        # compute the bits per pixel
        self.source._seek(self.ZERO)
        im = self.source.get_next_video_frame()
        im = num.fromstring(im.data,num.uint8)
        self.color_depth = len(im)/self.width/self.height
        if self.color_depth != 1 and self.color_depth != 3:
            raise ValueError( 'color_depth = %d, only know how to deal with color_depth = 1 or colr_depth = 3'%self.color_depth )
        self.bits_per_pixel = self.color_depth * 8

        # allocate the buffer
        self.buffer = num.zeros((self.height,self.width,self.buffersize),dtype=num.uint8)
        self.bufferts = num.zeros(self.buffersize)
        # put the first frame in it
        self.source._seek(self.ZERO)
        self.currframe = 0
        (frame,ts) = self.get_next_frame_and_reset_buffer()

        self._all_timestamps = None

    def get_all_timestamps(self):
        if self._all_timestamps is None:
            self._all_timestamps = []
            for f in range(self.n_frames):
                try:
                    (frame,ts) = self.get_frame(f)
                except:
                    continue
                self._all_timestamps.append(ts)
        return self._all_timestamps            

    def get_frame(self,framenumber):
        """Read frame from file and return as NumPy array."""

        if framenumber < 0: raise IndexError

        # have we already stored this frame?
        if framenumber >= self.bufferframe0 and framenumber < self.bufferframe1:
            off = num.mod(framenumber - self.bufferframe0 + self.bufferoff0,self.buffersize)
            if DEBUG: print "frame %d is in buffer at offset %d"%(framenumber,off)
            return (self.buffer[:,:,off].copy(),self.bufferts[off])

        # is framenumber the next frame to read in?
        if framenumber == self.currframe:
            if DEBUG: print "frame %d is the next frame, just calling get_next_frame"%framenumber
            return self.get_next_frame()

        # otherwise, we need to seek

        # find the last keyframe at or before framenumber
        lastkeyframe = num.floor(framenumber / self.keyframe_period)*self.keyframe_period
        lastkeyframe_s = lastkeyframe / self.fps + self.start_time

        if DEBUG: print "keyframe before frame %d is %d (ts = %f)"%(framenumber,lastkeyframe,lastkeyframe_s)

        didbeginbusy = False
        if params.interactive:
            wx.Yield()
            if not wx.IsBusy():
                didbeginbusy = True
                wx.BeginBusyCursor()

        # is the current frame before the desired frame and at or after the nearest keyframe?
        if framenumber > self.currframe and self.currframe >= lastkeyframe:

            if DEBUG: print "the best thing to do is just to read in frames until we get to desired %d (currframe = %d)"%(framenumber,self.currframe)

            # read into the buffer
            for f in range(self.currframe,framenumber+1):
                try:
                    (frame,ts) = self.get_next_frame()
                except:
                    if didbeginbusy:
                        wx.EndBusyCursor()
                    raise
            
            if didbeginbusy:
                wx.EndBusyCursor()

            return (frame,ts)
            
        # seek around til we find the right window

        if DEBUG: print "we need to seek"

        # make sure we don't move in opposite directions ever
        dirmoved = 0

        # convert from frame to time
        wantts = framenumber / self.fps + self.start_time
        if DEBUG: print "frame %d is equivalent to time stamp %f"%(framenumber,wantts)
        while True:
            # try seeking
            if DEBUG: print "start_time = %f"%self.start_time
            if DEBUG: print "trying to seek to %f (actually %f)"%(lastkeyframe_s,lastkeyframe_s-self.start_time)
            try: 
                self.source._seek(lastkeyframe_s-self.start_time)
            except:
                print 'seeking failed, aborting. tried to seek to %f (actually %f)'%(lastkeyframe_s,lastkeyframe_s-self.start_time)
            tscurr = self.source.get_next_video_timestamp()
            if tscurr is None:
                print 'seeking failed, aborting. timestamp returned is None'
                raise NoMoreFramesException

            if DEBUG: print "ended up seeking to %f"%tscurr

            # are we in the right ballpark?
            if wantts < tscurr:

                lastkeyframe_s = max(self.ZERO,lastkeyframe_s-self.keyframe_period_s)
                if DEBUG: print "still not back far enough, trying to seek to %f"%lastkeyframe_s
                dirmoved = -1
                continue

            elif wantts >= tscurr + self.keyframe_period_s:
                if dirmoved == -1:
                    # we went back and now we think we need to go forward
                    # so we're screwed up
                    # let's just stay
                    if DEBUG: print "uh oh -- we've already gone backward, and now we think we need to go forward, so we're just going to stay"
                    break
                lastkeyframe_s += self.keyframe_period_s
                dirmoved = 1
                continue
            else:
                # we're in the right ballpark
                if DEBUG: print "we're in the right ballpark"
                break

        # end loop searching for keyframe

        lastkeyframe = (tscurr - self.start_time)*self.fps
        self.currframe = int(lastkeyframe)

        (frame,ts) = self.get_next_frame_and_reset_buffer()

        # go forward until we get the right frame
        for f in range(self.currframe,framenumber+1):
            try:
                (frame,ts) = self.get_next_frame()
            except:
                if didbeginbusy:
                    wx.EndBusyCursor()
                raise

        if didbeginbusy:
            wx.EndBusyCursor()

        # return
        return (frame,ts)

    def get_next_frame_and_reset_buffer(self):

        # first frame stored in buffer
        self.bufferframe0 = self.currframe
        # frame after last frame stored in buffer
        self.bufferframe1 = self.currframe + 1
        # buffer will wrap around if frames are read in in-order
        # bufferoff0 is the location in the buffer of the first 
        # frame buffered
        self.bufferoff0 = 0
        (frame,ts) = self._get_next_frame_helper()
        self.buffer[:,:,0] = frame.copy()
        self.bufferts[0] = ts
        # bufferoff is the location in the buffer where we should 
        # write the next frame
        if self.buffersize > 1:
            self.bufferoff = 1
        else:
            self.bufferoff = 0

        # set the current frame
        self.currframe += 1
        self.prevts = ts

        return (frame,ts)

    def _get_next_frame_helper(self):
        ts = self.source.get_next_video_timestamp()
        if ts is None:
            if self.currframe - 1 < self.n_frames:
                self.n_frames = self.currframe - 1
            raise NoMoreFramesException
        im = self.source.get_next_video_frame()
        frame = num.fromstring(im.data,num.uint8)

        if self.color_depth == 1:
            frame.resize((im.height,im.width))
        else: # color_depth == 3
            frame.resize( (self.height, self.width*3) )
            tmp = frame.astype(float)
            tmp = tmp[:,2:self.width*3:3]*.3 + \
                tmp[:,1:self.width*3:3]*.59 + \
                tmp[:,0:self.width*3:3]*.11
            frame = tmp.astype(num.uint8)

        frame = num.flipud(frame)

        return (frame,ts)

    def get_next_frame(self):

        (frame,ts) = self._get_next_frame_helper()

        # store 

        # last frame stored will be 1 more
        self.bufferframe1 += 1

        # are we erasing the first frame?
        if self.bufferoff0 == self.bufferoff:
            self.bufferframe0 += 1
            if self.buffersize > 1:
                self.bufferoff0 += 1
            if DEBUG: print "erasing first frame, bufferframe0 is now %d, bufferoff0 is now %d"%(self.bufferframe0,self.bufferoff0)

        if DEBUG: print "buffer frames: [%d,%d), bufferoffset0 = %d"%(self.bufferframe0,self.bufferframe1,self.bufferoff0)
            
        self.buffer[:,:,self.bufferoff] = frame.copy()
        self.bufferts[self.bufferoff] = ts

        if DEBUG: print "read into buffer[%d], ts = %f"%(self.bufferoff,ts)

        self.bufferoff += 1

        # wrap around
        if self.bufferoff >= self.buffersize:
            self.bufferoff = 0

        if DEBUG: print "incremented bufferoff to %d"%self.bufferoff

        # remember current location in the movie
        self.currframe += 1
        self.prevts = ts

        if DEBUG: print "updated currframe to %d, prevts to %f"%(self.currframe,self.prevts)

        return (frame,ts)

    def _estimate_fps(self):

        if DEBUG: print 'Estimating fps'

        # seek to the start of the stream
        self.source._seek(self.ZERO)

        if DEBUG: print 'First seek succeeded'

        # initial time stamp
        ts0 = self.source.get_next_video_timestamp()
        ts1 = ts0

        if DEBUG: print 'initial time stamp = ' + str(ts0)

        # get the next frame and time stamp a bunch of times
        nsamples = 200
        if DEBUG: print 'nsamples = ' + str(nsamples)
        i = 0 # i is the number of frames we have successfully grabbed
        while True:
            im = self.source.get_next_video_frame()
            ts = self.source.get_next_video_timestamp()
            if DEBUG: print 'i = %d, ts = '%i + str(ts)
            if (ts is None) or num.isnan(ts) or (ts <= ts1):
                break
            i = i + 1
            ts1 = ts
            if i >= nsamples:
                break
        
        if ts1 <= ts0:
            raise ValueError( "Could not compute the fps in the compressed movie" )

        self.fps = float(i) / (ts1-ts0)
        if DEBUG: print 'Estimated frames-per-second = %f'%self.fps
        
    def _estimate_keyframe_period(self):

        if DEBUG: print 'Estimating keyframe period'

        self.source._seek(self.ZERO)

        ts0 = self.source.get_next_video_timestamp()

        if DEBUG: print 'After first seek, ts0 intialized to ' + str(ts0)

        i = 1 # i is the number of successful seeks
        foundfirst = False
        while True:

            # seek to the next frame
            self.source._seek(float(i)/self.fps)

            # get the current time stamp
            ts = self.source.get_next_video_timestamp()

            # did we seek past the end of the movie?
            if ts is None or num.isnan(ts):
                if foundfirst:
                    # we found one keyframe after start, use start
                    self.keyframe_period = i
                    self.keyframe_period_s = self.keyframe_period / self.fps
                else:
                    # then set keyframe period to be length of entire movie
                    self.keyframe_period = self.n_frames + 1
                    self.keyframe_period_s = self.duration_seconds + self.fps
                    if DEBUG: 'Only keyframe found at start of movie, setting keyframe_period = n_frames + 1 = %d, keyframe_period_s = duration_seconds + fps = %f'%(self.keyframe_period,self.keyframe_period_s)

                #raise ValueError( "Could not compute keyframe period in compressed video" )    
                return

            if ts > ts0:
                if foundfirst:
                    break
                else:
                    foundfirst = True
                    i0 = i
                    ts0 = ts

            i = i + 1

        if DEBUG: print "i = %d, i0 = %d"%(i,i0)
        self.keyframe_period = i - i0
        self.keyframe_period_s = self.keyframe_period / self.fps
        if DEBUG: print "Estimated keyframe period = " + str(self.keyframe_period)

    def get_n_frames( self ):
        return self.n_frames

    def get_width( self ):
        return self.width

    def get_height( self ):
        return self.height

    def seek(self,framenumber):

        ts = float(framenumber) / self.fps + self.start_time
        self.source._seek(ts)
        ts = self.source.get_next_video_timestamp()
        self.prevts = ts - 1. / self.fps
        self.currframe = (ts - self.start_time)*self.fps

        return self.currframe

def write_results_to_avi(movie,tracks,filename,f0=None,f1=None):

    nframes = len(tracks)
    if f0 is None:
        f0 = params.start_frame
    if f1 is None:
        f1 = nframes + params.start_frame - 1

    f0 -= params.start_frame
    f1 -= params.start_frame
    f0 = max(0,min(nframes-1,f0))
    f1 = max(0,min(nframes-1,f1))
    nframes_write = f1-f0+1

    # open the file for output
    outstream = open(filename,'wb')

    # write the header
    write_avi_header(movie,tracks,filename,outstream,f0,f1)

    # get the current location
    movilistloc = outstream.tell()

    # write the frames
    offsets = num.zeros(nframes_write)
    for i in range(f0,f1+1):
        if (i % 100) == 0:
            print 'Frame %d / %d'%(i,nframes_write)

        offsets[i-f0] = write_avi_frame(movie,tracks,i,outstream)

    # get offset relative to movilist
    offsets -= movilistloc + 4

    # write the index
    write_avi_index(movie,tracks,offsets,outstream,f0,f1)

    # close
    outstream.close()

def write_avi_index(movie,tracks,offsets,outstream,f0,f1):

    nframes = f1-f0+1
    idx1size = 8 + 16*nframes
    BYTESPERPIXEL = 3
    bytesperframe = movie.get_width()*movie.get_height()*BYTESPERPIXEL

    write_chunk_header('idx1',idx1size,outstream)

    for i in range(len(offsets)):
        outstream.write(struct.pack('4s','00db'))
        outstream.write(struct.pack('I',16))
        outstream.write(struct.pack('I',offsets[i]))
        outstream.write(struct.pack('I',bytesperframe))

def write_avi_frame(movie,tracks,i,outstream):

    height = movie.get_height()
    width = movie.get_width()
    BYTESPERPIXEL = 3
    bytesperframe = width*height*BYTESPERPIXEL

    if tracks is None:
        return
    if i >= len(tracks):
        return

    # global frame index
    j = params.start_frame + i

    # read in the video frame
    try:
        frame, last_timestamp = movie.get_frame(j)
    except (IndexError,NoMoreFramesException):
        return

    # get the current tracks
    ellipses = tracks[i]

    # get tails
    old_pts = []
    early_frame = int(max(0,i-params.tail_length))
    for dataframe in tracks[early_frame:i+1]:
        these_pts = []
        for ellipse in dataframe.itervalues():
            these_pts.append( (ellipse.center.x,ellipse.center.y,
                               ellipse.identity) )
        old_pts.append(these_pts)

    # draw on image
    bitmap,resize,img_size = annotate_bmp(frame,ellipses,old_pts,
                                            params.ellipse_thickness,
                                            [height,width])
    img = bitmap.ConvertToImage()
    # the image is flipped
    img = img.Mirror(True)
    img = img.GetData()

    # write chunktype
    outstream.write(struct.pack('4s','00db'))
    # write size of frame
    outstream.write(struct.pack('I',bytesperframe))

    # write frame
    offset = outstream.tell()
    outstream.write(img[::-1])
    pad = bytesperframe%2
    if pad == 1:
        outstream.write(struct.pack('B',0))
    return offset

def write_avi_header(movie,tracks,filename,outstream,f0,f1):

    # movie size
    BYTESPERPIXEL = 3
    nframes = f1-f0+1
    width = movie.get_width()
    height = movie.get_height()
    bytesperframe = width*height*BYTESPERPIXEL

    # chunk sizes if 0 frames
    avihsize = 64
    #strnsize = 8 + len(filename) + 1
    strllistsize = 116
    strhsize = 56
    strfsize = 48
    hdrllistsize = avihsize + strllistsize + 12
    movilistsize = 12
    idx1size = 8
    riffsize = hdrllistsize + movilistsize + idx1size
    # add in frames
    movilistsize += nframes * (4+4+bytesperframe+(bytesperframe%2))
    idx1size += nframes * (4*4)
    riffsize +=  nframes * (4+4+bytesperframe + 4*4 + (bytesperframe%2))
    ## add in strnsize
    #addon = strnsize + (strnsize%2)
    #riffsize += addon
    #hdrllistsize += addon
    #strllistsize += addon

    # write the RIFF chunk header
    write_chunk_header('RIFF',riffsize,outstream)
    # write AVI fourcc
    outstream.write(struct.pack('4s','AVI '))
    # write hdrl LIST
    write_list_header('hdrl',hdrllistsize-8,outstream)
    # write avih chunk
    write_chunk_header('avih',avihsize-8,outstream)

    ## write main avi header
    # microseconds per frame
    if hasattr(movie,'frame_delay_us'):
        microsecperframe = movie.frame_delay_us
    elif hasattr(movie.h_mov,'frame_delay_us'):
        microsecperframe = movie.h_mov.frame_delay_us
    else:
        microsecperframe = estimate_frame_delay_us(movie.h_mov)
    outstream.write(struct.pack('I',int(round(microsecperframe))))
    # maximum bytes per second
    framespersec = 1e6/microsecperframe
    bytespersec = framespersec*bytesperframe
    outstream.write(struct.pack('I',int(num.ceil(bytespersec))))
    # reserved
    outstream.write(struct.pack('I',0))
    # flags
    outstream.write(struct.pack('I',16))
    # number of frames
    outstream.write(struct.pack('I',nframes))
    # initial frame
    outstream.write(struct.pack('I',0))
    # number of streams
    outstream.write(struct.pack('I',1))
    # suggested buffer size
    outstream.write(struct.pack('I',bytesperframe))
    # width
    outstream.write(struct.pack('I',width))
    # height
    outstream.write(struct.pack('I',height))
    # frame rate
    outstream.write(struct.pack('2I',100,100*framespersec))
    # not sure -- start, length
    outstream.write(struct.pack('2I',0,0))

    # strl list
    write_list_header('strl',strllistsize-8,outstream)
    # strh chunk
    write_chunk_header('strh',strhsize-8,outstream)

    ## write stream header
    # FCC type
    outstream.write(struct.pack('4s','vids'))
    # FCC handler -- 'DIBS '
    outstream.write(struct.pack('I',0))
    # Flags
    outstream.write(struct.pack('I',0))
    # Reserved
    outstream.write(struct.pack('I',0))
    # Initial Frame
    outstream.write(struct.pack('I',0))
    # Frame rate
    outstream.write(struct.pack('2I',100,100*framespersec))
    # not sure -- start, length
    outstream.write(struct.pack('2I',0,0))
    # suggested buffer size
    outstream.write(struct.pack('I',bytesperframe))
    # quality
    outstream.write(struct.pack('I',7500))
    # not sure -- sample size
    outstream.write(struct.pack('I',0))

    # Write strf chunk
    write_chunk_header('strf',strfsize-8,outstream)

    ## Write bitmap header
    # Size
    outstream.write(struct.pack('I',40))
    # width
    outstream.write(struct.pack('I',width))
    # height
    outstream.write(struct.pack('I',height))
    # planes
    outstream.write(struct.pack('H',1))
    # bits per pixel
    outstream.write(struct.pack('H',24))
    # FourCC: DIBS
    outstream.write(struct.pack('I',0))
    # image size
    outstream.write(struct.pack('I',bytesperframe))
    # not sure
    outstream.write(struct.pack('4I',0,0,0,0))

    ## Write stream name chunk and data
    #write_chunk_header('strn',strnsize-8,outstream)
    #outstream.write(filename)
    #outstream.write(struct.pack('B',0))
    #if (len(filename)%2) == 1:
    #    outstream.write(struct.pack('B',0))

    # movi list
    write_list_header('movi',movilistsize,outstream)

def write_chunk_header(chunktype,chunksize,outstream):

    outstream.write(struct.pack('4sI',chunktype,chunksize))

def write_list_header(listtype,listsize,outstream):

    outstream.write(struct.pack('4sI4s','LIST',listsize,listtype))

def estimate_frame_delay_us(mov):

    if not hasattr(mov,'chunk_start'):
        return 0

    # go to beginning of first frame
    if mov.issbfmf:
        return .05*1e6
    else:
        mov.file.seek(mov.chunk_start)
        # read the first timestamp
        stamp0 = mov.get_next_timestamp()
        # go to the last frame
        mov.file.seek(mov.chunk_start+mov.bytes_per_chunk*(mov.n_frames-1))
        # read the last timestamp
        stamp1 = mov.get_next_timestamp()


        frame_delay_us = float(stamp1-stamp0)/float(mov.n_frames-1)*1e6
        return frame_delay_us
