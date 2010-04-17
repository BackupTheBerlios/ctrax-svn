import os, struct
import numpy as nx
import motmot.FlyMovieFormat.FlyMovieFormat as upstream_fmf

def load_movie(filename,check_integrity=False):
    """factory function to create StaticBackgroundFlyMovie or FlyMovie as appropriate"""
    if not isinstance(filename,basestring):
        # upstream_fmf.FlyMovie also supports file objects
        raise TypeError('filename must be a string')

    # get the extension
    tmp,ext = os.path.splitext(filename)
    if ext == '.sbfmf':
        return StaticBackgroundFlyMovie(filename,check_integrity=check_integrity)
    else:
        return upstream_fmf.FlyMovie(filename,check_integrity=check_integrity)

def FlyMovie(*args,**kw):
    """deprecated factory function that acts like old Ctrax.FlyMovie constructor"""
    import warnings
    warnings.warn("call Ctrax.load_movie() to load .fmf or .sbfmf movies", DeprecationWarning, stacklevel=2)
    return load_movie(*args,**kw)

class StaticBackgroundFlyMovie(upstream_fmf.FlyMovie):
    """an .sbfmf movie, loaded from a pre-existing file"""
    def __init__(self,file,check_integrity=False):
        # initialize base class
        upstream_fmf.FlyMovie.__init__(self,file,check_integrity=check_integrity)
        self.init_sbfmf()

    def init_sbfmf(self):
	#try:
        # read the version number
        format = '<I'
        nbytesver, = struct.unpack(format,self.file.read(struct.calcsize(format)))
        version = self.file.read(nbytesver)

        # read header parameters
        format = '<4IQ'
        nr,nc,self.n_frames,difference_mode,self.indexloc = \
	      struct.unpack(format,self.file.read(struct.calcsize(format)))

        # read the background image
        self.bgcenter = nx.fromstring(self.file.read(struct.calcsize('<d')*nr*nc),'<d')
        # read the std
        self.bgstd = nx.fromstring(self.file.read(struct.calcsize('<d')*nr*nc),'<d')

        # read the index
        ff = self.file.tell()
        self.file.seek(self.indexloc,0)
        self.framelocs = nx.fromstring(self.file.read(self.n_frames*8),'<Q')

	#except:
        #    raise InvalidMovieFileException('file could not be read')

	if version == "0.1":
	    self.format = 'MONO8'
	    self.bits_per_pixel = 8

	self.framesize = (nr,nc)
	self.bytes_per_chunk = None
	self.timestamp_len = struct.calcsize('<d')
        self.chunk_start = self.file.tell()
        self.next_frame = None
        self._all_timestamps = None # cache

    def _read_next_frame(self):
        format = '<Id'
        npixels,timestamp = struct.unpack(format,self.file.read(struct.calcsize(format)))
        idx = nx.fromstring(self.file.read(npixels*4),'<I')
        v = nx.fromstring(self.file.read(npixels*1),'<B')
        frame = self.bgcenter.copy()
        frame[idx] = v
        frame.shape = self.framesize
        return frame, timestamp

    def _read_next_timestamp(self):
        format = '<Id'
        self.npixelscurr,timestamp = struct.unpack(format,self.file.read(struct.calcsize(format)))
        return timestamp

    def get_all_timestamps(self):
        if self._all_timestamps is None:

            self._all_timestamps = []

            format = '<Id'
            l = struct.calcsize(format)
            for i in range(self.n_frames):
                self.seek(i)
                npixels,timestamp = struct.unpack(format,self.file.read(l))
                self._all_timestamps.append(timestamp)
            self.next_frame = None
            self._all_timestamps = nx.asarray(self._all_timestamps)
        return self._all_timestamps

    def seek(self,frame_number):
        if frame_number < 0:
            frame_number = self.n_frames + frame_number
        if frame_number < 0:
            raise IndexError("index out of range (n_frames = %d)"%self.n_frames)
        seek_to = self.framelocs[frame_number]
        self.file.seek(seek_to)
        self.next_frame = None
