# dump list of Matlab filenames into a text file so Matlab can read it
# JAB 5/12/11

import filelist

if __name__ == '__main__':

    files = filelist.CtraxFileList()

    # write original filenames
    fp = open( 'mat_name_dump.tmp', 'w' )

    for fname in files.mat_files():
        fp.write( '%s\n'%fname )

    fp.close()

    # write newly tracked filenames
    fp = open( 'mat_new_name_dump.tmp', 'w' )

    for fname in files.mat_tempfiles():
        fp.write( '%s\n'%fname )

    fp.close()
