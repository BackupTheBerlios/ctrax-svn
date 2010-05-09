"""houghcircles interface using ctypes.
"""
import numpy as num
from houghcircles_C import houghcircles_C

def houghcircles(x,y,w,binedgesa,bincentersb,bincentersr,hlib=None):

    # get data size
    npts = len(x)
    na = len(binedgesa)-1
    nb = len(bincentersb)
    nr = len(bincentersr)

    # allocate acc
    #acc = num.zeros(na*nb*nr,dtype=float)

    acc = houghcircles_C(x,y,w,binedgesa,bincentersb,bincentersr)

    #acc = acc.reshape([nr,nb,na])

    # since this code was written for matlab, order is backwards
    acc = acc.swapaxes(0,2)

    return acc
