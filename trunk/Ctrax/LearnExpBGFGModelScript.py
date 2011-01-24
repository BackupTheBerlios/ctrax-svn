import ExpBGFGModel
import re
import os
import os.path
import ParseExpDirs
import matplotlib.pyplot as plt

## parameters

mode = 'show'

# which experiments
protocol = 'CtraxTest20110111'
mindatestr = ''
maxdatestr = ''
linename = 'pB*'
rig = ''
plate = ''
bowl = ''
notstarted = False

# directory for learning data
resdir = '/groups/branson/bransonlab/tracking_data/olympiad/FlyBowl/CtraxTest20110111/LearnCtraxParams'

# files within the learning data directory
paramsFileStr = 'ExpBGFGModelParams.txt'
movieFileStr = 'movie.ufmf'
annFileStr = 'movie.ufmf.ann'
expdirsFileStr = 'expdirs_%s.txt'%protocol
outputFileStr = 'ExpBGFGModelResults.pickle'
matFileStr = 'ExpBGFGModelResults.mat'

# file to write experiment directory names to
expdirsFileName = os.path.join(resdir,expdirsFileStr)

# file containing parameters
paramsFileName = os.path.join(resdir,paramsFileStr)

# file to write results to
outputFileName = os.path.join(resdir,outputFileStr)

# mat file to write results to
matFileName = os.path.join(resdir,matFileStr)

# experiment directories corresponding to parameters
(expdirs,expdir_reads,expdir_writes,experiments) = \
    ParseExpDirs.getExpDirs(protocol=protocol,
                            mindatestr=mindatestr,
                            maxdatestr=maxdatestr,
                            linename=linename,
                            rig=rig,
                            plate=plate,
                            bowl=bowl,
                            notstarted=notstarted)

if mode == 'learn':

    print "expdirs = " + str(expdirs)

    fid = open(expdirsFileName,"w")

    for expdir in expdir_reads:
        fid.write('%s\n'%expdir)

    fid.close()

    print 'execute the following command:'
    print 'python ExpBGFGModel.py' + \
        ' -f ' + expdirsFileName + \
        ' -p ' + paramsFileName + \
        ' -m ' + movieFileStr + \
        ' -a ' + annFileStr + \
        ' -o ' + outputFileName + \
        ' --mat ' + matFileName

elif mode == 'show':
    
    model = ExpBGFGModel.ExpBGFGModel(picklefile=outputFileName)
    model.show()
    for i in range(len(expdirs)):
        moviename = os.path.join(expdir_reads[i],movieFileStr)
        model.showtest(moviename=moviename)
        savename = os.path.join(resdir,'ExpBGFGModel_SampleFrames_%s.png'%expdirs[i])
        plt.savefig(savename,format='png')
