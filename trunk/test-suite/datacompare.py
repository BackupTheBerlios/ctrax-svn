# statistics for Ctrax performance comparisons
# JAB 5/14/11

import os
import time

from scipy.io import loadmat, savemat

class CtraxDataComparator:
    def __init__( self, old_matfiles, new_matfiles ):
        """Read data from old and new MAT-files."""

        self.olddata = []
        for oldfile in old_matfiles:
            data = loadmat( oldfile )
            self.olddata.append( {'x':data['x'],
                                  'y':data['y'],
                                  'theta':data['theta'],
                                  'a':data['a'],
                                  'b':data['b']} )

        self.newdata = []
        for newfile in new_matfiles:
            data = loadmat( newfile )
            self.newdata.append( {'x':data['x'],
                                  'y':data['y'],
                                  'theta':data['theta'],
                                  'a':data['a'],
                                  'b':data['b']} )


    def save_data( self, filelist, runtimes ):
        """Save formatted comparisons to disk."""

        stat_filenames = filelist.stat_files()

        for olddata, newdata, fname, runtime in zip( self.olddata, self.newdata,
                                                     stat_filenames, runtimes ):
            savemat( fname, {'truedata':olddata, 'newdata':newdata, 'runtime':runtime} )

        # write stats filenames (which are timestamped!)
        time_str = time.strftime( '%Y-%m-%d_%H-%M-%S' )
        dump_name = os.path.join( filelist.tmp_dir, 'stats_' + time_str + '.tmp' )

        fp = open( dump_name, 'w' )
        for fname in stat_filenames:
            fp.write( '%s\n'%fname )
        fp.close()

        return dump_name
