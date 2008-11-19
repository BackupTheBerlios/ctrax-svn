from medmadpyx import medmadpyx
import numpy as num

def medmad(filename,nframes,nr,nc,nframesskip,bytesperchunk,headersize):

    print "headersize = %d"%headersize

    (medl,madl) = medmadpyx(filename,nframes,nr,nc,nframesskip,
                            bytesperchunk,headersize)

    med = num.array(medl,copy=True)
    mad = num.array(madl,copy=True)

    print "mad range = [%f,%f]"%(num.min(mad),num.max(mad))
    
    med.shape = (nr,nc)
    mad.shape = (nr,nc)
    
    return (med,mad)
