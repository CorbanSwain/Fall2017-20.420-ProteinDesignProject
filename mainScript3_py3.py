from datetime import datetime
from pyrosetta import *
from pyrosetta.rosetta.protocols.simple_moves import *
from pyrosetta.rosetta.protocols.relax import FastRelax
import os

def madeGlobal(varName):
    return varName in globals()


def initialize():
    if not madeGlobal('didInit'):
        init()
        global didInit, defaultScorefxn
        didInit = True
        defaultScorefxn = get_fa_scorefxn()

initialize()
    
def now(format=0):
    if format == 0:
        return datetime.now().strftime('[%H:%M:%S %m-%d-%y]')
    else:
        return datetime.now().strftime('%y%m%d')

def dprint(text):
    text = ' {}  '.format(text)
    print('[{}] {}'.format(now(),text.center(70,'*')))

def printScore(pose,title,scorefxn=defaultScorefxn):
    title = title + ' Score'
    print('{} --> {:9.5f}'.format(title.lalign(30), scorefxn(pose_mzn)))
    
def loadInPose(fileName):
    dprint('Loading In `{}` and Creating Pose'.format(filename))
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

def smallNShearMove(pose,\
                    scorefxn=defaultScorefxn,\
                    repeats=50,\
                    n_moves=5,\
                    kT=1.0,\
                    angle=0):
    dprint('Running Small and Shear Movers {:2d} Times'.format(repeats))
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

def mkDir(directory):
    if not os.path.isdir(directory):
        os.makedirs(directory)

def isFile(fileName):
    return os.path.isfile(fileName)

def setupCaches():
    if not madeGlobal('cacheDir'):
        global cacheDir
        cacheDir = 'AlgorithmCache'
    mkDir(cacheDir)

    if not madeGlobal('pdbCache'):
        global pdbCache
        pdbCache = os.path.join(cacheDir,'PDBs')
    mkdir(pdbCache)

    global sesCache
    if madeGlobal('sesCache'):
        if now(1) in sesCache:
            return
    sesCache = '{}_session_log.txt'.format{now(1)}
    sesCache = os.path.join(cacheDir,sesCache)
    if not isFile(sesCache):
        with  open(sesCache,'w') as f:
            f.writelines('{}\n\n'.format(sesCache))
            f.writelines('Protein Design Algorithm Trials, Session Log\nBegins: {}\n\n'.format(now())) 

def main():
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
    n1 = 5
    angle = 25
    
    for i in range(n1):
        dprint('Beginning Loop # {:2d}/{:2d}'.format(i+1,n1))
        angle = angle - 5
        smallNShearMove(minPose,angle=angle)

    dprint('Done With Loop')
    pymol.apply(minPose)

    printScore(startPose,'Original')
    printScore(fRelaxPose,'Fast')
    printScore(minPose,'Minimization')

main()

 



