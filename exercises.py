from datetime import datetime
from pyrosetta import *
from pyrosetta.rosetta.protocols.simple_moves import *
from pyrosetta.rosetta.protocols.relax import FastRelax

init()

def now():
    return datetime.now().strftime('[%H:%M:%S %m-%d-%y]')

def dprint(text):
    text = '  {}  '.format(text)
    print('\n\n{}\n'.format(text.center(100,'*')))

dprint('Loading In PDB, Creating Pose')
complexFile = 'new_3vi8_complex.pdb'
params = ['LG.params']
pose = Pose()
generate_nonstandard_residue_set(pose, params)
pose_from_file(pose,complexFile)
relaxed_pose = Pose()
relaxed_pose.assign(pose)

print('Num Residues: {:d}'.format(pose.total_residue()))

dprint('Loading Score Function')
scorefxn = get_fa_scorefxn()

dprint('[{}] Running Fast Relax'.format(now()))
fast_relax = FastRelax()
fast_relax.set_scorefxn(scorefxn) 
fast_relax.apply(relaxed_pose)
dprint('[{}] Done Fast Relax'.format(now()))

pose_mzn = Pose()
pose_mzn.assign(pose)

dprint('Setting up Move Map')
kT = 1.0
n_moves = 5
movemap = MoveMap()
movemap.set_bb_true_range(130, 185)

dprint('Setting Up Min Mover')
min_mv = MinMover()
min_mv.movemap(movemap)
min_mv.score_function(scorefxn)

(n_1, n_2) = 100, 5

def repMoverWithAngle(angle):
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
    mc = MonteCarlo(pose_mzn, scorefxn, kT)

    dprint('Building Sequence Mover')
    seq_mv = SequenceMover()
    seq_mv.add_mover(small_mv)
    seq_mv.add_mover(min_mv)
    seq_mv.add_mover(shear_mv)
    seq_mv.add_mover(min_mv)

    dprint('Building Trial and Repeat Movers')
    trial_mv = TrialMover(seq_mv, mc)
    rep_mv = RepeatMover(trial_mv, n_1)
    return rep_mv

angle = 25    

dprint('[{}] Beginning For Loop'.format(now()))
for i in range(n_2):
    print('[{}] Beginning Loop # {:d}'.format(now(),i+1))
    repMoverWithAngle(angle).apply(pose_mzn)
    angle = angle - 5

dprint('[{}] Done With Loop'.format(now()))

print('Original Score    : {:05f}'.format(scorefxn(pose)))
print('Fast Relax Score  : {:05f}'.format(scorefxn(relaxed_pose)))
print('Minimization Score: {:05f}'.format(scorefxn(pose_mzn)))

dprint('Sending Poses to PyMOL')
pymol = PyMOLMover()
pose.pdb_info().name('original')
relaxed_pose.pdb_info().name('relaxed')
pose_mzn.pdb_info().name('minimized')

pymol.apply(pose)
pymol.apply(relaxed_pose)
pymol.apply(pose_mzn)
 
