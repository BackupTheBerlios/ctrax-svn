import os, glob

for filename in glob.glob(os.path.join('..','*.py')):
    fidin = open(filename,'r')
    (head,tail) = os.path.split(filename)
    tmpfilename = os.path.join(head,'rpr_tmp_' + tail)
    fidout = open(tmpfilename,'w')
    lineno = 0
    ischange = False
    print "\n\nProcessing file = %s:"%filename
    for line in fidin:
        lineno+=1
        is_pkg_resources = line.find('pkg_resources')
        # no pkg_resources, then direct copy
        if is_pkg_resources == -1:
            fidout.write(line)
        else:
            is_import = line.find('import ')
            if is_import == -1:
                is_resource_filename = line.find('resource_filename')
                if is_resource_filename == -1:
                    print '**[%s:%d] Could not fix:\n  %s'%(filename,lineno,line)
                    fidout.write(line)
                else:
                    quotestart = max(line.find('"',is_pkg_resources),line.find("'",is_pkg_resources))
                    quoteend = max(line.find('"',quotestart+1),
                                   line.find("'",quotestart+1))
                    parenend = line.find(')',is_pkg_resources)
                    if quotestart == -1 or quoteend == -1 or parenend == -1:
                        print '**[%s:%d]: Found resource_filename, but could not fix:\n  %s'%(filename,lineno,line)
                    else:
                        fidout.write(line[:is_pkg_resources] + line[quotestart:quoteend+1] + line[parenend+1:])
                        ischange = True
                        print '[%s:%d] Removed resource_filename command:\n  %s'%(filename,lineno,line)

            else:
                tmp = line.strip()
                if tmp[0] == '#':
                    print '[%s:%d] Import is already commented:\n  %s'%(filename,lineno,line)
                    fidout.write(line)
                else:
                    fidout.write('# ' + line)
                    ischange = True
                    print '[%s:%d] Commented:\n  %s'%(filename,lineno,line)
    fidin.close()
    fidout.close()
    if ischange:
        print "Replacing old %s"%filename
        os.remove(filename)
        os.rename(tmpfilename,filename)
    else:
        print "No changes made to %s"%filename
        os.remove(tmpfilename)
