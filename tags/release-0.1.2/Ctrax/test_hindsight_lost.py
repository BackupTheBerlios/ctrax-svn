import ellipsesk as ell
import numpy as num
import matchidentities
import hindsight
from params import params
import estconncomps as est
import pylab as mpl

colors = ['r','g','b','y']
tracks = []

params.dampen = 0.

# frame 0
tracks.append(ell.TargetList())
# just one ellipse with id = 0
tracks[0][0] = ell.Ellipse(10.,10.,1.,2.,num.pi/3,0.,0)
tracks[0][0].compute_area()
# frame 1
# just one ellipse with id = 0
tracks.append(ell.TargetList())
tracks[1][0] = tracks[0][0].copy()
tracks[1][0].x += 10.
tracks[1][0].y += 10.
# frame 2
# detection missed at pred
# new ellipse at another location
tracks.append(ell.TargetList())
#tracks[2][1] = ell.Ellipse(20.,20.,1.,2.,num.pi/3,0.,1)
prev = tracks[0].copy()
curr = tracks[1].copy()
pred = matchidentities.cvpred(prev,curr)
# frame 3
# no ellipses, but predicted to be at pred
tracks.append(ell.TargetList())
prev = curr
curr = pred
pred = matchidentities.cvpred(prev,curr)
# frame 4
# new ellipse at pred
tracks.append(ell.TargetList())
prev = curr
curr = pred
pred = matchidentities.cvpred(prev,curr)
pred0 = pred[0].copy()
pred0.identity = 1
tracks[4][1] = pred0
# frame 5
# new ellipse at pred
tracks.append(ell.TargetList())
prev = curr
curr = pred
pred = matchidentities.cvpred(prev,curr)
pred0 = pred[0].copy()
pred0.identity = 1
tracks[5][1] = pred0
# frame 6
# new ellipse at pred
tracks.append(ell.TargetList())
prev = curr
curr = pred
pred = matchidentities.cvpred(prev,curr)
pred0 = pred[0].copy()
pred0.identity = 1
tracks[6][1] = pred0

params.nids = 2

bw = []
bounds = [0,100,0,100]
nr = bounds[1]-bounds[0]
nc = bounds[3]-bounds[2]
for track in tracks:
    tmp = num.zeros((nr,nc),dtype=bool)
    for ellipse in track.itervalues():
        tmp |= est.ellipsepixels(ellipse,bounds)
    bw.append(tmp)

dfore = []
for b in bw:
    dfore.append(b.astype(float))

class BG:
    def __init__(self,dfore,bw):
        self.dfore = dfore
        self.bw = bw
    def sub_bg(self,frame):
        return (self.dfore[frame],self.bw[frame])

#for t in range(len(tracks)):
#    mpl.imshow(bw[t])
#    for (id,e) in tracks[t].iteritems():
#        est.drawellipse(e,colors[id])
#    mpl.show()

bg = BG(dfore,bw)
    
# try to fix
hind = hindsight.Hindsight(tracks,bg)
hind.initialize_milestones()
print 'milestones: '
print hind.milestones

print 'tracks: '
print tracks

for t in range(2,len(tracks)):
    hind.fixerrors(t)

for t in range(len(tracks)):

    mpl.imshow(bw[t])
    for [id,e] in tracks[t].iteritems():
        est.drawellipse(e,colors[id])
    mpl.show()


