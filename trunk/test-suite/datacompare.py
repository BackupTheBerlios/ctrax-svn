# statistics for Ctrax performance comparisons
# JAB 5/14/11

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


    def save_data( self, stat_filenames, runtimes ):
        """Save formatted comparisons to disk."""
        pass


    
