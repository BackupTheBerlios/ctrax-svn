from params import params
params.interactive = False
import ExpBGFGModel
import re
import os
import os.path
import ParseExpDirs
import sys
import matplotlib.pyplot as plt

## parameters

mode = 'show'

# which experiments
protocol = '20110407'
expdirs_in = '/groups/branson/bransonlab/projects/olympiad/FlyBowlCtrax/20110407/LearnCtraxParams/expdirs.txt'
mindatestr = ''
maxdatestr = ''
linename = '.*'
rig = ''
plate = ''
bowl = ''
notstarted = False

print "*** PROTOCOL %s ***\n"%protocol

# directory for learning data
resdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlCtrax/%s/LearnCtraxParams'%protocol

# directory containing parameters
paramsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlCtrax/%s'%protocol

# files within the learning data directory
paramsFileStr = 'ExpBGFGModelParams.txt'
movieFileStr = 'movie.ufmf'
annFileStr = 'movie.ufmf.ann'
expdirsFileStr = 'expdirs.txt'
outputFileStr = 'ExpBGFGModelResults.pickle'
matFileStr = 'ExpBGFGModelResults.mat'


# file to write experiment directory names to
expdirsFileName = os.path.join(resdir,expdirsFileStr)

# file containing parameters
paramsFileName = os.path.join(paramsdir,paramsFileStr)

# file to write results to
outputFileName = os.path.join(resdir,outputFileStr)

# mat file to write results to
matFileName = os.path.join(resdir,matFileStr)

# experiment directories corresponding to parameters
if expdirs_in == "":
    (expdirs,expdir_reads,expdir_writes,experiments) = \
        ParseExpDirs.getExpDirs(protocol=protocol,
                                mindatestr=mindatestr,
                                maxdatestr=maxdatestr,
                                linename=linename,
                                rig=rig,
                                plate=plate,
                                bowl=bowl,
                                notstarted=notstarted,
                                subreadfiles=[movieFileStr,annFileStr])
else:
    fid = open(expdirs_in,"r")
    expdirs = []
    experiments = []
    for expdir in fid:
        expdir = expdir.strip()
        if expdir == "":
            continue
        expdirs.append(expdir)
        [exp,success] = ParseExpDirs.parseExpDir(expdir)
        experiments.append(exp)
    expdir_reads = expdirs
    expdir_writes = expdirs
    fid.close()

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
        print moviename
        model.showtest(moviename=moviename)
        savename = os.path.join(resdir,'ExpBGFGModel_SampleFrames_%s.png'%os.path.basename(expdirs[i]))
        plt.savefig(savename,format='png')

elif mode == 'change':
    model = ExpBGFGModel.ExpBGFGModel(picklefile=outputFileName)
