#!/usr/bin python3
# mainScript3_py3.py

# from datetime import datetime
from pyrosetta import *
from datetime import datetime
from pyrosetta.rosetta.protocols.simple_moves import *
from pyrosetta.rosetta.protocols.relax import FastRelax
import os

def madeGlobal(varName):
    return varName in globals()

def now(formatNum=0):
    if formatNum == 0:
        return datetime.now().strftime('%H:%M:%S %m-%d-%y')
    elif formatNum == 1:
        return datetime.now().strftime('%y%m%d')
    else:
        return datetime.now().strftime('%H:%M:%S')

def mkDir(directory):
    if not os.path.isdir(directory):
        os.makedirs(directory)

def isFile(fileName):
    return os.path.isfile(fileName)

def initialize():
    if not madeGlobal('didInit'):
        init()
        global didInit, defaultScorefxn
        didInit = True
        defaultScorefxn = get_fa_scorefxn()
initialize()

def setupCaches():
    if not madeGlobal('cacheDir'):
        global cacheDir
        cacheDir = 'AlgorithmCache'
    mkDir(cacheDir)

    if not madeGlobal('pdbCache'):
        global pdbCache
        pdbCache = os.path.join(cacheDir,'PDBs')
    mkDir(pdbCache)

    global sesCache
    if madeGlobal('sesCache'):
        if now(1) in sesCache:
            return
    sesCache = '{}_session_log.txt'.format(now(1))
    sesCache = os.path.join(cacheDir,sesCache)
    if not isFile(sesCache):
        with open(sesCache,'w') as f:
            f.writelines('{}\n\n'.format(sesCache))
            f.writelines('Protein Design Algorithm Trials, Session Log\nBegins: {}\n\n'\
                         .format(now())) 
setupCaches()
    
def log(text=''):
    with open(sesCache,'a') as f:
        f.writelines('[{}] >  {}\n'.format(now(2),text))

def logBegin():
    with open(sesCache,'a') as f:
        f.writelines('\n\n[{}] >  {}\n'.format(now(2),' Beginning Script '.center(70,'v')))

def logEnd():
    with open(sesCache,'a') as f:
        f.writelines('[{}] >  {}\n\n'.format(now(2),' Script Complete! '.center(70,'^')))
        
def dprint(text):
    log(text)
    text = ' {} '.format(text)
    print('[{}] {}'.format(now(),text.center(70,'*')))
        
def printScore(pose,title,scorefxn=defaultScorefxn):
    title = title + ' Score'
    output = '{} --> {:11.5f}'.format(title.ljust(30), scorefxn(pose))
    log(output)
    print(output)
    
def loadInPose(fileName):
    dprint('Loading In `{}` and Creating Pose'.format(fileName))
    params = ['LG.params']
    pose = Pose()
    generate_nonstandard_residue_set(pose, params)
    pose_from_file(pose,fileName)
    return pose

def poseFrom(pose):
    newPose = Pose()
    newPose.assign(pose)
    return newPose

def namePose(pose,nameStr):
     pose.pdb_info().name(nameStr)

def fastRelax(pose,scorefxn=defaultScorefxn):
    dprint('Beginning Fast Relax')
    fast_relax = FastRelax()
    fast_relax.set_scorefxn(scorefxn) 
    fast_relax.apply(pose)

def minMove(pose,\
            scorefxn=defaultScorefxn,\
            repeats=10,
            kT=5):
    dprint('Perfoming A Minimization Loop')
    log('Num Repeats: {:4d} | kT: {:9.2f}'.\
        format(repeats,kT))
    movemap = MoveMap()
    movemap.set_bb(True)
    
    min_mv = MinMover()
    min_mv.movemap(movemap)
    min_mv.score_function(scorefxn)

    mc = MonteCarlo(pose,scorefxn,kT)

    trial_mv = TrialMover(min_mv,mc)
    rep_mv = RepeatMover(trial_mv,repeats)

    rep_mv.apply(pose)
    
    
    
def smallNShearMove(pose,\
                    scorefxn=defaultScorefxn,\
                    repeats=50,\
                    n_moves=5,\
                    kT=1.0,\
                    angle=0):
    dprint('Running Small and Shear Movers'.format(repeats))
    log('Num Repeats: {:4d} | Num Sml/Shr Moves: {:4d} | kT: {:9.2f} | Angl Rng: {:4.1f}'.\
        format(repeats,n_moves,kT,angle))
    movemap = MoveMap()
    movemap.set_bb(True)

    min_mv = MinMover()
    min_mv.movemap(movemap)
    min_mv.score_function(scorefxn)

    small_mv = SmallMover(movemap, kT, n_moves)
    shear_mv = ShearMover(movemap, kT, n_moves)
    if angle:
        small_mv.angle_max('H', angle)
        small_mv.angle_max('E', angle)
        small_mv.angle_max('L', angle)
        shear_mv.angle_max('H', angle)
        shear_mv.angle_max('E', angle)
        shear_mv.angle_max('L', angle)

    mc = MonteCarlo(pose, scorefxn, kT)

    seq_mv = SequenceMover()
    seq_mv.add_mover(small_mv)
    seq_mv.add_mover(min_mv)
    seq_mv.add_mover(shear_mv)
    seq_mv.add_mover(min_mv)

    trial_mv = TrialMover(seq_mv, mc)
    rep_mv = RepeatMover(trial_mv, repeats)

    rep_mv.apply(pose)

def createPyMolMover():
    pymol = PyMOLMover()
    return pymol

def main():
    logBegin()
    pymol = createPyMolMover()

    startPose = loadInPose('new_3vi8_complex.pdb')
    print('Num Residues: {:d}'.format(startPose.total_residue()))
    pymol.apply(startPose)
    namePose(startPose,'original')

    fRelaxFile = '3vi8_complex_fastRelaxed.pdb'
    fRelaxFile = os.path.join(pdbCache,fRelaxFile)
    if os.path.isfile(fRelaxFile):
        fRelaxPose = loadInPose(fRelaxFile)
    else:
        fRelaxPose = poseFrom(startPose)
        fastRelax(fRelaxPose)
        fRelaxPose.dump_pdb(fRelaxFile)
        namePose(fRelaxPose,'orig_relaxed')
    pymol.apply(fRelaxPose)
    
    minPose = poseFrom(startPose)

    angles = [20.0, 18.0, 18.0, 10.0] + [4.0,] * 4
    kT1s = [1000.0, 800.0, 800.0, 500.0] + [20,] * 4
    annealTime = len(angles)
    kT2 = 10.0
    repeats = 10
    n_moves = 5

    dprint('Beginning Small Shear Anneal Loop')
    n1 = len(kT1) * 4
    log('Cycles: {:2d} | kT: {:3.1f}'.format(n1,kT2))
    for i in range(n1):
        dprint('Beginning Loop # {:2d}/{:2d}'.format(i+1,n1))
        mc = MonteCarlo(minPose, defaultScorefxn, kT2)
        ind = i % annealTime
        angle = angles[ind]
        kT = kT1s[ind]
        smallNShearMove(minPose,\
                        angle=angle, kT=kT1[i % len(kT1)],\
                        repeats=repeats, n_moves=n_moves)
        mc.boltzmann(minPose)
        printScore(minPose,'Iteration # {:2d}'.format(i+1))
        
    dprint('Done With Loop')

    n2 = 200
    minMove(minPose,repeats=n2)
    printScore(minPose,'After {} Min Cycles'.format(n2))
    pymol.apply(minPose)

    log()
    printScore(startPose,'Original')
    printScore(fRelaxPose,'Fast')
    printScore(minPose,'Minimization')
    logEnd()
    
main()

 



