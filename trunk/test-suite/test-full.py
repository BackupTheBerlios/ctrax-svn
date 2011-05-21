#!/usr/bin/env python

# main script for full Ctrax tracking tests
# JAB 5/10/11

import os
import time

import datacompare
import filelist

def main():
    # get filenames
    files = filelist.CtraxFileList()

    # run Ctrax on each file
    tmp_files = files.movie_tempfiles()
    runtimes = []
    for filename in tmp_files:
        starttime = time.time()

        cmd = "Ctrax-script.py --Interactive=False --Input=%s"%filename
        #os.system( cmd )

        runtimes.append( time.time() - starttime )

    # get Matlab to make readable versions of the test suite's original (fixed) MAT-files
    os.system( "matlab -nodisplay -nojvm -nosplash -r 'make_test_data; exit'" )

    # compare original MAT-data with newly generated MAT-data
    comparator = datacompare.CtraxDataComparator( files.resaved_mat_files(), files.resaved_mat_tempfiles() )
    out_filename = comparator.save_data( files, runtimes )

    # run visualization in Matlab
    print "wrote tracking data to ", out_filename
    #os.system( "matlab -r \"plot_test_data_comparisons( '" + out_filename + "' )\"" )
    

if __name__ == '__main__':
    main()
