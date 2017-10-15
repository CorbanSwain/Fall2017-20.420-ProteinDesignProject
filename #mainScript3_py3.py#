#!/usr/bin python3
# mainScript3_py3.py

# from datetime import datetime
from pyrosetta import *
from datetime import datetime
from pyrosetta.toolbox import generate_resfile_from_pose
from pyrosetta.rosetta.protocols.simple_moves import *
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.protocols.moves import AddPyMOLObserver
import os
import re

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
        global didInit, defaultScorefxn, proteinName
        proteinName = '3vi8_BPDE'
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

def makeSequenceGlobals():
    if madeGlobal('madeSeqGlobals'):
        if madeSeqGlobals == True:
            return
    global pocketResNums, firstResNum, oneIndexedRes,\
        allAa, aaGroupings, adnlAa, conservMutMap, liberalMutMap, madeSeqGlobals
    madeSeqGlobals = False
    pocketResNums = [241, 247, 250, 251, 254, 255, 272, 273, 275, 276, 279, 280,\
                     321, 325, 330, 332, 333, 334, 339, 343, 344, 354, 355, 358]
    firstResNum = 202
    oneIndexedRes = [i - firstResNum + 1 for i in pocketResNums]
    allAa = ['A','C','D','E','F','G','H','I','K','L',\
             'M','N','P','Q','R','S','T','V','W','Y']
    print(len(allAa))
    aaGroupings = [['A','I','L','M','V'],\
                   ['F','W','Y'],\
                   ['N','C','Q','T','S'],\
                   ['D','E'],\
                   ['R','H','K'],\
                   ['G'],\
                   ['P']]
    adnlAa = [['F','W','Y','G'],\
              ['A','I','L','M','V','G'],\
              ['Y','W','G'],\
              ['R','H','K','G'],\
              ['D','E','G'],\
              ['A','S'],\
              ['G']]
    conservMutMap = {}
    liberalMutMap = {}
    for i,aaList in enumerate(aaGroupings):
        for j,aa in enumerate(aaList):
            temp_list = list(aaList)
            del temp_list[j]
            conservMutMap[aa] = temp_list
            liberalMutMap[aa] = temp_list + adnlAa[i]
    madeSeqGlobals = True
makeSequenceGlobals()

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
            kT=1):
    dprint('Perfoming A Minimization Loop')
    log(' `--> Num Repeats: {:4d} | kT: {:9.2f}'.\
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
                    repeats=50,\
                    n_moves=5,\
                    kT=1.0,\
                    angle=0):
    dprint('Running Small and Shear Movers'.format(repeats))
    log(' `--> Num Repeats: {:4d} | Num Sml/Shr Moves: {:4d} | kT: {:9.2f} | Angl Rng: {:4.1f}'.\
        format(repeats,n_moves,kT,angle))
    movemap = MoveMap()
    movemap.set_bb(True)

    min_mv = MinMover()
    min_mv.movemap(movemap)
    min_mv.score_function(defaultScorefxn)

    small_mv = SmallMover(movemap, kT, n_moves)
    shear_mv = ShearMover(movemap, kT, n_moves)
    if angle:
        small_mv.angle_max('E', angle)
        small_mv.angle_max('L', angle)
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

def smallNShearAnnealLoop(pose,cycles,kT1):
    angles = [5, 3, 3, 2] + [0.5,] * 4
    kT2s = [100.0, 50.0, 50.0, 30.0] + [5,] * 4
    annealTime = len(angles)
    seqRepeats = 5
    numSmallShearRepeats = 5
    n1 = len(kT2s) * annealCycles
    
    dprint('Beginning Small/Shear Anneal Loop')

    log(' `--> Cycles: {:2d} | kT: {:3.1f}'.format(n1,kT2))
    for i in range(n1):
        dprint('Beginning Loop # {:2d}/{:2d}'.format(i+1,n1))
        mc = MonteCarlo(pose, defaultScorefxn, kT1)
        ind = i % annealTime
        angle = angles[ind]
        kT2 = kT2s[ind]
        smallNShearMove(pose,\
                        angle=angle, kT=kT2,\
                        repeats=seqRepeats, n_moves=numSmallShearRepeats)
        mc.boltzmann(pose)
        printScore(pose,'Iteration # {:2d}'.format(i+1))
    
def createPyMolMover():
    pymol = PyMOLMover()
    return pymol

def generateResfile(pose,name=proteinName):
    fext = 'resfile'
    fname = '.'.join([name,fext])
    generate_resfile_from_pose(pose, fname)

def mutateResfile(pose,\
                  nameout='MUT',\
                  namein=proteinName,\
                  residues=pocketResNums,\
                  liberal=False):
    fext = 'resfile'
    with open('.'.join([namein, fext])) as rfile:
        text = rfile.read()

    if liberal: mutDict = liberalMutMap
    else: mutDict = conservMutMap
    reg = r'(?P<num>{})(?P<chain>\s+\w\s+)(?P<cmds>\w+)'.\
          format('|'.join(['{:3d}'.format(r) for r in residues]))
    def subFun(m):
        num = int(m.group('num').strip())
        chain = m.group('chain').strip()
        poseNum = pose.pdb_info().pdb2pose(chain,num)
        aa = pose.residue(poseNum).name1()
        mutRes = mutDict[aa]
        return '{:3d}  {}  PIKAA {}'.format(num,chain,''.join(mutRes))
         
    new_text = re.sub(reg,subFun,text)
    fnameOut = '.'.join(['{}_{}'.format(namein, nameout), fext]) 
    with open(fnameOut, 'w+') as new_rfile:
        new_rfile.write(new_text)
    return fnameOut

def main():
    logBegin()
    pymol = createPyMolMover()

    startPose = loadInPose('new_3vi8_complex.pdb')
    print('Num Residues: {:d}'.format(startPose.total_residue()))
    namePose(startPose,'original')
    pymol.apply(startPose)

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
    
    minPose = poseFrom(fRelaxPose)
    namePose(minPose,'minimized')

    annealCycles = 2
    kT = 10.0
    smallNShearAnnealLoop(minPose,annealCycles,kT)

    n2 = 200
    kT3 = 0.1
    minMove(minPose,repeats=n2,kT=kT3)
    printScore(minPose,'After {} Min Cycles'.format(n2))
    pymol.apply(minPose)

    printScore(startPose,'Original')
    printScore(fRelaxPose,'Fast Relax')
    printScore(minPose,'Minimization Protocol')
    logEnd()
    
main()

 



