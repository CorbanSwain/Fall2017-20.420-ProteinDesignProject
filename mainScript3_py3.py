from datetime import datetime
from pyrosetta import *
from pyrosetta.rosetta.protocols.simple_moves import *
from pyrosetta.rosetta.protocols.relax import FastRelax

def initialize():
    if 'didInit' in globals():
        return
    init()
    global didInit, defaultScorefxn
    didInit = True
    defaultScorefxn = get_fa_scorefxn()

initialize()
    
def now():
    return datetime.now().strftime('[%H:%M:%S %m-%d-%y]')

def dprint(text):
    text = ' {}  '.format(text)
    print('[{}] {}'.format(now(),text.center(70,'*')))

def loadInPose():
    dprint('Loading In PDB, Creating Pose')
    complexFile = 'new_3vi8_complex.pdb'
    params = ['LG.params']
    pose = Pose()
    generate_nonstandard_residue_set(pose, params)
    pose_from_file(pose,complexFile)
    return pose

def poseFromPose(pose):
    newPose = Pose()
    newPose.assign(pose)
    return newPose

def namePose(pose,nameStr):
     pose.pdb_info().name(nameStr)

def fastRelax(pose,scorefxn=defaultScorefxn):
    dprint('[{}] Beginning Fast Relax'.format(now()))
    fast_relax = FastRelax()
    fast_relax.set_scorefxn(scorefxn) 
    fast_relax.apply(pose)
    dprint('[{}] Done Fast Relax'.format(now()))

def smallNShearMove(pose,\
                    scorefxn=get_fa_scorefxn(),\
                    repeats=50,\
                    n_moves=5,\
                    kT=1.0,\
                    angle=15):
    movemap = MoveMap()
    movemap.set_bb(True)
    min_mv = MinMover()
    min_mv.movemap(movemap)
    min_mv.score_function(scorefxn)
    dprint('Setting Up Simple Movers')
    small_mv = SmallMover(movemap, kT, n_moves)
    shear_mv = ShearMover(movemap, kT, n_moves)
    small_mv.angle_max('H', angle)
    small_mv.angle_max('E', angle)
    small_mv.angle_max('L', angle)
    shear_mv.angle_max('H', angle)
    shear_mv.angle_max('E', angle)
    shear_mv.angle_max('L', angle)
    dprint('Setting Up MonteCarlo')
    mc = MonteCarlo(p, scorefxn, kT)
    dprint('Building Sequence Mover')
    seq_mv = SequenceMover()
    seq_mv.add_mover(small_mv)
    seq_mv.add_mover(min_mv)
    seq_mv.add_mover(shear_mv)
    seq_mv.add_mover(min_mv)
    dprint('Building Trial and Repeat Movers')
    trial_mv = TrialMover(seq_mv, mc)
    rep_mv = RepeatMover(trial_mv, repeats)
    rep_mv.apply(pose)
    return rep_mv   

def main():
    pymol = PyMOLMover()

    startPose = loadInPose()
    print('Num Residues: {:d}'.format(startPose.total_residue()))
    pymol.apply(startPose)
    namePose(startPose,'original')

    fRelaxPose = poseFromPose(startPose)
    namePose(fRelaxPose,'orig_relaxed')
    fastRelax(fRelaxPose)

    dprint('Beginning For Loop')
    n1 = 5
    for i in range(n_1):
        dprint('Beginning Loop # {:02d}'.format(i+1))
        smallNShearMove(minPose)
        angle = angle - 5

    dprint('[{}] Done With Loop'.format(now()))
    pymol.apply(pose_mzn)

    print('Original Score    : {:05f}'.format(scorefxn(pose)))
    print('Fast Relax Score  : {:05f}'.format(scorefxn(relaxed_pose)))
    print('Minimization Score: {:05f}'.format(scorefxn(pose_mzn)))





 



