import wx
from wx import xrc
import motmot.wxvalidatedtext.wxvalidatedtext as wxvt
import numpy as num
from params import params

import pkg_resources # part of setuptools
RSRC_FILE = pkg_resources.resource_filename( __name__, "chooseorientations.xrc" )

class ChooseOrientations:

    def __init__(self,parent,targets,interactive=True):

        self.interactive = interactive
        self.targets = targets

        if self.interactive:
            rsrc = xrc.XmlResource( RSRC_FILE )
            self.frame = rsrc.LoadFrame(parent,"orientationframe")

        self.InitControlHandles()
        self.InitializeValues()
        self.BindCallbacks()


    def InitControlHandles(self):

        if self.interactive:
            self.weight_input = self.control('weight')
            self.max_weight_input = self.control('max_weight')

    def InitializeValues(self):

        self.weight = params.velocity_angle_weight
        self.max_weight = params.max_velocity_angle_weight
        if self.weight is None:
            self.weight = .5/params.max_jump
        if self.max_weight is None:
            self.max_weight = 1
        if self.interactive:
            self.weight_input.SetValue('%.3f'%self.weight)
            self.max_weight_input.SetValue('%.3f'%self.max_weight)

    def BindCallbacks(self):

        if self.interactive:
            wxvt.setup_validated_float_callback(self.weight_input,
                                                xrc.XRCID('weight'),
                                                self.ValidateWeight,
                                                pending_color=params.wxvt_bg)
            wxvt.setup_validated_float_callback(self.max_weight_input,
                                                xrc.XRCID('weight'),
                                                self.ValidateWeight,
                                                pending_color=params.wxvt_bg)

    def control(self,ctrlname):
        return xrc.XRCCTRL(self.frame,ctrlname)

    def ValidateWeight(self,evt):
        self.weight = float(self.weight_input.GetValue())
        self.max_weight = float(self.max_weight_input.GetValue())
        if self.weight < 0:
            self.weight = 0
        if self.max_weight < 0:
            self.max_weight = 0
        if self.max_weight > 1:
            self.max_weight = 1
        
        self.weight_input.SetValue('%.3f'%self.weight)
        self.max_weight_input.SetValue('%.3f'%self.max_weight)

    def anglemod(self,theta):

        return ( ( (theta+num.pi)%(2.*num.pi) )-num.pi )

    def angledist(self,theta1,theta2):
        
        return abs( ( (theta1-theta2+num.pi)%(2.*num.pi) ) - num.pi)

    def ChooseOrientations(self):

        params.velocity_angle_weight = self.weight
        params.max_velocity_angle_weight = self.max_weight

        N = len(self.targets)
        startframes = num.zeros(params.nids,dtype=int)
        startframes[:] = -1
        endframes = num.zeros(params.nids,dtype=int)
        endframes[:] = -1
        keystostart = set(range(params.nids))
        keystoend = set([])
        allkeys = set(range(params.nids))
        for t in range(N):
            keys = set(self.targets[t].keys())
            newstarts = keystostart & keys
            for i in newstarts:
                startframes[i] = t
            keystostart -= newstarts
            keystoend = keystoend | newstarts
            keys = allkeys - keys
            newends = keystoend & keys
            for i in newends:
                endframes[i] = t
            keystoend -= newends
        for i in keystoend:
            endframes[i] = N

        for i in range(params.nids):
            if startframes[i] < endframes[i]-1:
                self.ChooseOrientationsPerID(i,startframes[i],endframes[i])

    def ChooseOrientationsPerID(self,id,startframe,endframe):
        
        # we will set the angle to theta_t = phi_t + s_t * pi
        # we want to choose s_t to minimize
        # \sum_t cost(s_t|s_{t-1})
        # cost(s_t|s_{t-1}) = [(1 - w(||v_t||^2))*d(\theta_t,\theta_{t-1}) +
        #                      w(||v_t||^2)*d(\theta_t,angle(v_t))]
        # where w(||v_t||^2) = \min{1, c*||v_t||^2}
        # we will find the most likely states s_t using the recursion
        # cost_t(s_t) = min_{s_{t-1}} { cost_{t-1}(s_{t-1}) + cost(s_t|s_{t-1})
        
        N = endframe - startframe
        
        # allocate space for storing the optimal path
        stateprev = num.zeros((N-1,2),dtype=bool)
        
        # allocate space for computing costs
        tmpcost = num.zeros(2)
        costprevnew = num.zeros(2)
        
        # initialize first frame
        costprev = num.zeros(2)
        
        # compute iteratively
        for tloc in range(1,N):

            t = tloc + startframe
	    
	    # compute velocity
            xcurr = self.targets[t][id].center.x
            ycurr = self.targets[t][id].center.y
	    xprev = self.targets[t-1][id].center.x
	    yprev = self.targets[t-1][id].center.y
	    vx = xcurr - xprev
	    vy = ycurr - yprev
	    # compute angle of velocity
	    velocityangle = num.arctan2(vy,vx)
            # compute weight for velocity term
	    w = num.minimum(float(params.max_velocity_angle_weight),params.velocity_angle_weight*num.sqrt(vx**2+vy**2))
            # compute for both possible states
            for scurr in [0,1]:
                
                thetacurr = self.targets[t][id].angle + scurr*num.pi
                
                # try both previous states
                for sprev in [0,1]:
                    
                    thetaprev = self.targets[t-1][id].angle + sprev*num.pi
                    
                    costcurr = (1.-w)*self.angledist(thetaprev,thetacurr) + \
                               w*self.angledist(thetacurr,velocityangle)
                    
                    tmpcost[sprev] = costprev[sprev] + costcurr

                # end loop over sprev
                
                # choose the minimum
                sprev = num.argmin(tmpcost)
                
                # set pointer for path
                stateprev[tloc-1,scurr] = sprev
                
                # set cost
                costprevnew[scurr] = tmpcost[sprev]

            # end loop over scurr

            costprev[:] = costprevnew[:]
	
	# end loop over frames

	# choose the best last state
	scurr = num.argmin(costprev)
	if scurr == 1:
	    self.targets[endframe-1][id].angle += num.pi
	    self.targets[endframe-1][id].angle = self.anglemod(self.targets[endframe-1][id].angle)

        # choose the best states
        for tloc in range(N-2,-1,-1):
            
            t = tloc + startframe
            
            scurr = stateprev[tloc,scurr]
            if scurr == 1:
                self.targets[t][id].angle += num.pi
                self.targets[t][id].angle = self.anglemod(self.targets[t][id].angle)

