class Avi:
    """Read uncompressed AVI movies."""
    def __init__( self, filename ):
        print "avi"
        self.file = open( filename, 'r' )
        self.read_header()

    def __del__( self ):
        """Existence possibly unnecessary, if Python does this cleanup for us."""
        if hasattr( self, 'file' ):
            self.file.close()

    def read_header( self ):
        # RIFF header
        file_type, riff_size = struct.unpack( '4sI', self.file.read( 8 ) )
        print "file_type = " + str(file_type)
        print "riff_size = " + str(riff_size)
        assert file_type == 'RIFF'
        stream_type = self.file.read( 4 )
        print "stream_type = " + str(streamp_type)
        assert stream_type == 'AVI '
        header_list, header_listsize, header_listtype = \
                     struct.unpack( '4sI4s', self.file.read( 12 ) )
        print "header_list = " + str(header_list)
        print "header_listsize = " + str(header_listsize)
        assert header_list == 'LIST' and header_listtype == 'hdrl'
        #size 4588 (fmf 1222)
        avi_str, avi_note = struct.unpack( '4sI', self.file.read( 8 ) )
        assert avi_str == 'avih'

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
        if fccType != 'vids':
            raise TypeError( "stream type must be video" )
        if fccHandler != 'DIB ':
            raise TypeError( "video must be uncompressed" )
        if stream_buf_size != 0:
            if hasattr( self, 'buf_size' ):
                assert self.buf_size == stream_buf_size
            else:
                self.buf_size = stream_buf_size

        # bitmap header
        bmp_str, bmp_note = struct.unpack( '4sI', self.file.read( 8 ) )
        assert bmp_str == 'strf'

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

        # skip extra header crap
        movie_list = ''
        movie_listtype = ''
        while movie_list != 'LIST':
            s = ''
            EOF_flag = False
            while s.find( 'movi' ) < 0 and not EOF_flag:
                p = self.file.tell()
                s = self.file.read( 128 )
                if s == '': EOF_flag = True
            if EOF_flag: break
            self.file.seek( p )
            self.file.read( s.find( 'movi' ) - 8 )
            movie_list, movie_listsize, movie_listtype = \
                        struct.unpack( '4sI4s', self.file.read( 12 ) )
        assert movie_list == 'LIST' and movie_listtype == 'movi'

        # beginning of data blocks
        self.data_start = self.file.tell()

    def get_frame( self, framenumber ):
        """Read frame from file and return as NumPy array."""
        if framenumber < 0: raise IndexError
        if framenumber >= self.n_frames: raise NoMoreFramesException

        # read frame from file
        self.file.seek( self.data_start + (self.buf_size+8)*framenumber )
        this_frame_id, frame_size = struct.unpack( '4sI', self.file.read( 8 ) )
        if frame_size != self.buf_size:
            raise ValueError( "movie must be uncompressed" )
        if not hasattr( self, 'frame_id' ):
            self.frame_id = this_frame_id
        elif this_frame_id != self.frame_id:
            raise ValueError( "error seeking frame start: unknown data header" )
        frame_data = self.file.read( frame_size )

        # make frame into numpy array
        frame = num.fromstring( frame_data, num.uint8 )
        if frame.size == self.width*self.height:
            frame.resize( (self.height, self.width) )
        elif frame.size == self.width*self.height*3:
            #frame.resize( (self.height, self.width, 3) )
            #frame = imops.to_mono8( 'RGB8', frame )
            raise TypeError( "movie must be grayscale" )
        else:
            raise ValueError( "frame size %d doesn't make sense: movie must be 8-bit grayscale"%(frame.size) )

        # make up a timestamp based on the file's stated framerate
        if self.frame_delay_us != 0:
            stamp = framenumber * self.frame_delay_us / 1e6
        elif self.time_scale != 0:
            stamp = framenumber * self.data_rate / float(self.time_scale)
        else:
            stamp = framenumber / 24
            # should raise warning or error here?
            # this might screw up playback, at least

        return frame, stamp

    def get_n_frames( self ): return self.n_frames
    def get_width( self ): return self.width
    def get_height( self ): return self.height
