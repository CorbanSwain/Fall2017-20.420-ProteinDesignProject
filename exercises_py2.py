from datetime import datetime
from pyrosetta import *
from pyrosetta.rosetta.protocols.simple_moves import *
from pyrosetta.rosetta.protocols.relax import FastRelax

init()

def now():
    return datetime.now().strftime("%I:%M%:S%p on %B %d, %Y")

def dprt(text):
    print '%s  %s  %s' % ('*'*30,text,'*'*30)

dprt('Loading In PDB, Creating Pose')
complexFile = 'new_3vi8_complex.pdb'
params = ['LG.params']
pose = Pose()
generate_nonstandard_residue_set(pose, params)
pose_from_file(pose,complexFile)
relaxed_pose = Pose()
relaxed_pose.assign(pose)

print 'Num Residues: %d' % (pose.total_residue())

dprt('Loading Score Function')
scorefxn = get_fa_scorefxn()

dprt('[%s] Running Fast Relax' % (now()))
fast_relax = FastRelax()
fast_relax.set_scorefxn(scorefxn) 
# fast_relax.apply(relaxed_pose)
dprt('[%s] Done Fast Relax' % (now()))

pose_mzn = Pose()
pose_mzn.assign(pose)

dprt('Setting up Move Map')
movemap = MoveMap()
movemap.set_bb(True) 

dprt('Setting Up Min Mover')
min_mv = MinMover()
min_mv.movemap(movemap)
min_mv.score_function(scorefxn) 

(n_1, n_2) = 100, 5 
kT = 1.0
n_moves = 5 

dprt('Setting Up Simple Movers')
small_mv = SmallMover(movemap, kT, n_moves)
shear_mv = ShearMover(movemap, kT, n_moves)

def moverWithAngle(angle):
    small_mv.angle_max('H', angle)
    small_mv.angle_max('E', angle)
    small_mv.angle_max('L', angle)
    shear_mv.angle_max('H', angle)
    shear_mv.angle_max('E', angle)
    shear_mv.angle_max('L', angle)
    
    dprt('Setting Up MonteCarlo')
    mc = MonteCarlo(pose_mzn, scorefxn, kT)

    dprt('Building Sequence Mover')
    seq_mv = SequenceMover()
    seq_mv.add_mover(small_mv)
    seq_mv.add_mover(min_mv)
    seq_mv.add_mover(shear_mv)
    seq_mv.add_mover(min_mv)

    dprt('Building Trial and Repeat Movers')
    trial_mv = TrialMover(seq_mv, mc)
    rep_mv = RepeatMover(trial_mv, n_1)
    return rep_mv

angle = 25
dprt('[%s] Beginning For Loop' % (now()))

for i in range(n_2):
    print '[%s] Beginning Loop # %d' % (now(),i+1)
    # moverWithAngle(angle).apply(pose_mzn)
    angle = angle - 5
    print 'Energy: %05f' % scorefxn(pose_mzn)
    
dprt('[%s] Done With Loop' % (now()))

print 'Original Score    : %05f' % (scorefxn(pose))
print 'Fast Relax Score  : %05f' % (scorefxn(relaxed_pose))
print 'Minimization Score: %05f' % (scorefxn(pose_mzn))

dprt('Sending Poses to PyMOL')
pymol = PyMOLMover()
pose.pdb_info().name('original')
relaxed_pose.pdb_info().name('relaxed')
pose_mzn.pdb_info().name('minimized')

pymol.apply(pose)
pymol.apply(relaxed_pose)
pymol.apply(pose_mzn)
 
