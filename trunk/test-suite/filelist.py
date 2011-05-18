# list of files for Ctrax tracking testing
# JAB 5/10/11

import os
import shutil
import stat
import tempfile
import time


class CtraxFileList:
    def __init__( self ):
        """Make a list of filenames to be used in testing."""
        
        use_datatypes = ['roach_fmf']#, 'sample_sbfmf1', 'sample_sbfmf2',
#                         'flybowl_ufmf1', 'flybowl_ufmf2', 'flybowl_ufmf3', 'flybowl_ufmf4',
#                         'flybowl_ufmf5', 'flybowl_ufmf6', 'flybowl_ufmf7', 'flybowl_ufmf8']
        
        self.base_dir = os.path.join( os.environ['HOME'], 'data', 'ctrax-data', 'test-suite' )

        self.tmp_dir = os.path.join( tempfile.gettempdir(), 'Ctrax-test-data' )
        if not os.access( self.tmp_dir, os.F_OK ):
            os.mkdir( self.tmp_dir )

        roach_fmf_dir = 'tietz'
        roach_fmf = ['movie20100707_141914_trimmed.fmf',
                     'movie20100707_142727_trimmed.fmf',
                     'movie20100707_143623_trimmed.fmf',
                     'movie20100707_144439_trimmed.fmf',
                     'movie20100707_145259_trimmed.fmf',
                     'movie20100707_150215_trimmed.fmf',
                     'movie20100707_151155_trimmed.fmf',
                     'movie20100707_152028_trimmed.fmf',
                     'movie20100707_152835_trimmed.fmf',
                     'movie20100707_153630_trimmed.fmf',
                     'movie20100707_154448_trimmed.fmf']


        sample_sbfmf1_dir = os.path.join( 'samples-distributed',
                                          'example-sbfmfvideo-ann-mat-5M5F5min-20071009_163231' )
        sample_sbfmf1 = ['movie20071009_163231.sbfmf']
        
        sample_sbfmf2_dir = os.path.join( 'samples-distributed',
                                          'example-sbfmfvideo-ann-mat-25M25F5min-20071009_164618' )
        sample_sbfmf2 = ['movie20071009_164618.sbfmf']


        flybowl_ufmf1_dir = os.path.join( 'Fixed_AR', 'EP00005_rc3',
                                          'pBDPGAL4U_TrpA_Rig1Plate10BowlA_20110304T154712' )
        flybowl_ufmf1 = ['movie.ufmf']

        flybowl_ufmf2_dir = os.path.join( 'Fixed_AR', 'EP00005_rc3',
                                          'pBDPGAL4U_TrpA_Rig1Plate10BowlB_20110302T155958' )
        flybowl_ufmf2 = ['movie.ufmf']
        
        flybowl_ufmf3_dir = os.path.join( 'Fixed_AR', 'EP00005_rc3',
                                          'pBDPGAL4U_TrpA_Rig1Plate10BowlC_20110303T111924' )
        flybowl_ufmf3 = ['movie.ufmf']


        flybowl_ufmf4_dir = os.path.join( 'Fixed_AR', 'EP00005_rc8',
                                          'GMR_64F08_AE_01_TrpA_Rig2Plate14BowlC_20110310T112011' )
        flybowl_ufmf4 = ['movie.ufmf']

        flybowl_ufmf5_dir = os.path.join( 'Fixed_AR', 'EP00005_rc8',
                                          'GMR_70F07_AE_01_TrpA_Rig1Plate10BowlC_20110310T130507' )
        flybowl_ufmf5 = ['movie.ufmf']

        flybowl_ufmf6_dir = os.path.join( 'Fixed_AR', 'EP00005_rc8',
                                          'GMR_74H07_AE_01_TrpA_Rig2Plate14BowlB_20110316T152501' )
        flybowl_ufmf6 = ['movie.ufmf']

        flybowl_ufmf7_dir = os.path.join( 'Fixed_AR', 'EP00005_rc8',
                                          'pBDPGAL4U_TrpA_Rig1Plate10BowlB_20110310T151713' )
        flybowl_ufmf7 = ['movie.ufmf']

        flybowl_ufmf8_dir = os.path.join( 'Fixed_AR', 'EP00005_rc8',
                                          'pBDPGAL4U_TrpA_Rig2Plate14BowlC_20110316T155136' )
        flybowl_ufmf8 = ['movie.ufmf']


        # make master file list
        self.files = []
        self.barefiles = []
        for dtype in use_datatypes:
            dirname = os.path.join( self.base_dir, eval( dtype + '_dir' ) )
            
            filenames = eval( dtype )
            for fname in filenames:
                self.files.append( os.path.join( dirname, fname ) )
                self.barefiles.append( fname )

        # test file existence
        for fname, annname, matname in zip( self.files,
                                            self.ann_files(),
                                            self.mat_files() ):
            if not os.access( fname, os.R_OK ):
                print "**%s not readable, skipping"%fname
                self.files.remove( fname )
            elif not os.access( annname, os.R_OK ):
                print "**annotation file %s not readable, skipping"%annname
                self.files.remove( fname )
            elif not os.access( matname, os.R_OK ):
                print "**MAT-file %s not readable, skipping"%matname
                self.files.remove( fname )
            

    def movie_files( self ): return self.files


    def ann_files( self ):
        """Make an annotation filename for each movie file."""
        
        ann_filelist = []
        for fname in self.files:
            ann_filelist.append( fname + '.ann' )
        return ann_filelist


    def mat_files( self, filelist=None ):
        """Make a Matlab filename for each movie file."""

        if filelist is None:
            filelist = self.files
        mat_filelist = []

        for fname in filelist:
            if fname.endswith( '.fmf' ):
                shortname = fname[:fname.index( '.fmf' )]
            elif fname.endswith( '.sbfmf' ):
                shortname = fname[:fname.index( '.sbfmf' )]
            elif fname.endswith( '.ufmf' ):
                shortname = fname[:fname.index( '.ufmf' )]
            elif fname.endswith( '.avi' ):
                shortname = fname[:fname.index( '.avi' )]
            else:
                shortname = fname

            # look for fixed versions of these
            matname = shortname + '.mat'
            head, tail = os.path.split( matname )
            fixed_tail = 'fixed_' + tail
            fixed_name = os.path.join( head, fixed_tail )
                
            if os.access( fixed_name, os.R_OK ):
                mat_filelist.append( fixed_name )
            else:
                mat_filelist.append( matname )

        return mat_filelist


    def resaved_mat_files( self, mat_files=None ):
        """Make filenames for MAT-files which have been resaved by Matlab."""

        if mat_files is None:
            mat_files = self.mat_files()
        resaved_files = []
        for fname in mat_files:
            resaved_files.append( fname[:fname.index( '.mat' )] + '_forpython.mat' )

        return resaved_files
    

    def make_movie_tempfiles( self ):
        """Copy movie files to temporary location for tracking."""

        self.tempfiles = []

        for fname, barename in zip( self.files, self.barefiles ):
            # look in temp location for file as new or newer, else copy
            tempname = os.path.join( self.tmp_dir, barename )
            if os.access( tempname, os.R_OK ):
                orig_stat = os.stat( fname )
                tmp_stat = os.stat( tempname )
                if tmp_stat[stat.ST_MTIME] > orig_stat[stat.ST_MTIME]:
                    print "%s has recent modification -- recopying to temp dir"%barename
                    os.remove( tempname )
                    shutil.copy2( fname, tempname )
            else:
                print "%s not found -- copying to temp dir"%barename
                shutil.copy2( fname, tempname )

            self.tempfiles.append( tempname )
    

    def movie_tempfiles( self ):
        """Get a list of temporary files for each movie file."""
        if not hasattr( self, 'tempfiles' ):
            self.make_movie_tempfiles()

        return self.tempfiles


    def mat_tempfiles( self ):
        """Get a list of MAT-file names for each temporary movie file."""

        return self.mat_files( self.movie_tempfiles() )


    def resaved_mat_tempfiles( self ):
        """Get a list of resaved MAT-file names for each temporary movie file."""

        return self.resaved_mat_files( self.mat_tempfiles() )


    def stat_files( self ):
        """Make a list of MAT-file names to save comparison statistics in."""

        if not hasattr( self, 'stat_names' ):
            filenames = self.mat_files()
            self.stat_names = []
            for fname in filenames:
                time_str = time.strftime( '%Y-%m-%d_%H-%M-%S' )
                self.stat_names.append( fname[:fname.index( '.mat' )] + '_' + time_str + '_comparison.mat' )

        return self.stat_names
