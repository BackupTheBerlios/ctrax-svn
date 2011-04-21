import re
import os
import os.path

def parseExpDir(expdir):

    expr = '(?P<line>[^/]+)_(?P<effector>[^_]*)_Rig(?P<rig>[0-9]+)Plate(?P<plate>[0-9]+)Bowl(?P<bowl>[A-Z]+)_(?P<notstarted>(notstarted_)?)(?P<date>[^_/]+)/?$'
    m = re.search(expr,expdir)
    d = dict()
    if m == None:
        success = False
        d['line'] = ''
        d['effector'] = ''
        d['rig'] = ''
        d['plate'] = ''
        d['bowl'] = ''
        d['notstarted'] = False
        d['date'] = ''
    else:
        success = True
        d['line'] = m.group('line')
        d['effector'] = m.group('effector')
        d['rig'] = m.group('rig')
        d['plate'] = m.group('plate')
        d['bowl'] = m.group('bowl')
        d['notstarted'] = m.group('notstarted') == 'notstarted'
        d['date'] = m.group('date')

    return (d,success)

def getExpDirs(protocol = '',
               mindatestr = '',
               maxdatestr = '',
               linename = '',
               rig = '',
               plate = '',
               bowl = '',
               notstarted = False,
               subreadfiles = [],
               subwritefiles = []):

    if protocol == 'scratched_polycarbonate_CtraxTest20101118':
        rootreaddir = '/groups/branson/bransonlab/tracking_data/olympiad/FlyBowl/polycarbonate_scratched'
        rootwritedir = '/groups/branson/home/bransonk/tracking/data/olympiad/FlyBowl/CtraxTest20101118'

    elif protocol == '20101118':
        rootreaddir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data/'
        rootwritedir = '/groups/branson/home/bransonk/tracking/data/olympiad/FlyBowl/CtraxTest20101118'

    elif protocol == '20110111':
        rootreaddir = '/groups/branson/bransonlab/tracking_data/olympiad/FlyBowl/CtraxTest20110111'
        rootwritedir = '/groups/branson/bransonlab/tracking_data/olympiad/FlyBowl/CtraxTest20110111'

    elif protocol == '20110202':
        rootreaddir = '/groups/branson/bransonlab/tracking_data/olympiad/FlyBowl/CtraxTest20110202'
        rootwritedir = '/groups/branson/bransonlab/tracking_data/olympiad/FlyBowl/CtraxTest20110202'
    else:
        rootreaddir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data/'
        rootwritedir = rootreaddir


    expdir_reads = set([f for f in os.listdir(rootreaddir) if os.path.isdir(os.path.join(rootreaddir, f))])

    expdir_writes = set([f for f in os.listdir(rootwritedir) if os.path.isdir(os.path.join(rootreaddir, f))])

    expdirs1 = expdir_reads.intersection(expdir_writes)

    expdir_reads = []
    expdir_writes = []
    expdirs = []
    experiments = []

    for expdir in expdirs1:

        [parsedcurr,success] = parseExpDir(expdir)
        if success == False:
            continue

        # check date range
        if mindatestr != '' and mindatestr > parsedcurr['date']:
            continue

        if maxdatestr != '' and maxdatestr < parsedcurr['date']:
            continue

        # check line name
        if linename != '' and re.search(linename,parsedcurr['line']) is None:
            continue

        # check rig
        if rig != '' and re.search(rig,parsedcurr['rig']) is None:
            continue

        # check plate
        if plate != '' and re.search(plate,parsedcurr['plate']) is None:
            continue

        # check bowl
        if bowl != '' and re.search(bowl,parsedcurr['bowl']) is None:
            continue

        # check notstarted
        if notstarted is not None and notstarted != parsedcurr['notstarted']:
            continue

        # check subfiles
        doexist = True
        for f in subreadfiles:
            if not os.path.exists(os.path.join(rootreaddir,expdir,f)):
                doexist = False
                break
        if not doexist:
            continue
        for f in subwritefiles:
            if not os.path.exists(os.path.join(rootwritedir,expdir,f)):
                doexist = False
                break
        if not doexist:
            continue

        expdir_reads.append(os.path.join(rootreaddir,expdir))
        expdir_writes.append(os.path.join(rootwritedir,expdir))
        expdirs.append(expdir)
        experiments.append(parsedcurr)

    return (expdirs,expdir_reads,expdir_writes,experiments)
