# ProdeinDesignScripts

from pyrosetta import *
import datetime

def now():
    return datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y")

outFile = open('output.txt','a+')
outFile.write('\n\n' + '*'*30 + '\n')
outFile.write(now())
outFile.write('\n' + '*'*30 + '\n')

def debugHeader(text):
    print('\n{0} {1} {0}\n'.format('*'*15,text))     

def initialize():
    init()
    
def load():
    import subprocess
    ligandFile = 'pah_aligned'
    debugHeader('CALLING OPEN BABEL')
    subprocess.call('babel {0}.pdb  {0}.mdl'.format(ligandFile), shell=True)
    debugHeader('MOLFILE 2 PARAMS')
    subprocess.call('molfile2params2/molfile_to_params.py --clobber {}.mdl'.format(ligandFile),\
                    shell=True)

def importPdb():
    debugHeader('Running import.Pdb')
    global pose, scorefxn, complexFile, pose2
    complexFile = 'new_3vi8_complex'
    params_list = ['LG.params']
    pose = Pose()
    debugHeader('Generating Nonstandard Res Set')
    generate_nonstandard_residue_set(pose, params_list)
    debugHeader('Generating Pose from File')
    pose_from_file(pose,'{}.pdb'.format(complexFile))
    print(pose)
    scorefxn = get_fa_scorefxn()

def sendtopymol():
    global pymol
    pymol = pyrosetta.PyMOLMover()
    pymol.apply(pose)

    print('my score: %s' % scorefxn(pose))
    outFile.write(str(scorefxn(pose)) + '\n')
    pymol.send_energy(pose)

def showhbonds():
    pymol.send_hbonds(pose)

def findmimima():
    min_mover = rosetta.protocols.simple_moves.MinMover()
    move_map = MoveMap()
    move_map.set_bb_true_range(40,130)
    min_mover.movemap(move_map)
    min_mover.score_function(scorefxn)

    pose_min_move = Pose()
    pose_min_move.assign(pose)

    observer = rosetta.protocols.moves.AddPyMOLObserver(pose, True)

    min_mover.apply(pose_min_move)

    print('original energy: %.2f' % scorefxn(pose))
    
    print('energy after min_mover: %.2f' % scorefxn(pose_min_move))

    outFile.write('original energy: %.2f\n' % scorefxn(pose))
    outFile.write('energy after min_mover: %.2f\n' % scorefxn(pose_min_move))
    
def fastrelax():
    fast_relax =rosetta.protocols.relax.FastRelax()
    fast_relax.set_scorefxn(scorefxn)

    pose_fast_relax = Pose()
    pose_fast_relax.assign(pose)

    pose_fast_relax.pdb_info().name('fast_relax')
    print(pose.pdb_info().name)
    print(pose_fast_relax.pdb_info().name())

    observer =rosetta.protocols.moves.AddPyMOLObserver(pose_fast_relax, True)

    fast_relax.apply(pose_fast_relax)

    print('finished!')

def fastrelax2():

    fast_relax  =rosetta.protocols.relax.FastRelax()
    fast_relax.set_scorefxn(scorefxn)

    pose_fast_relax = Pose()
    pose_fast_relax.assign(pose)
    outFile.write('Beginning Fast Relax: {}\n'.format(now()))
    fast_relax.apply(pose_fast_relax)

    print('original energy: %.2f' % scorefxn(pose))
    print('energy after fast relax: %.2f' % scorefxn(pose_fast_relax))

    outFile.write('Ending Fast Relax: {}\n'.format(now()))
    outFile.write('original energy: %.2f\n' % scorefxn(pose))
    outFile.write('energy after fast relax: %.2f\n' % scorefxn(pose_fast_relax))

    pose_fast_relax.pdb_info().name('fast_relax_final')
    pymol.apply(pose_fast_relax)
    pymol.send_energy(pose_fast_relax)


def mutataRes():
    pose_pack_mover = Pose()
    pose_pack_mover.assign(pose)
    pose_pack_mover.pdb_info().name('pack_mover')

    task = standard_packer_task(pose_pack_mover)
    task.temporarily_fix_everything()
    task.temporarily_set_pack_residue(124, True)
    task.temporarily_set_pack_residue(123, True)

    pack_mover = rosetta.protocols.simple_moves.PackRotamersMover(scorefxn, task)

    observer = rosetta.protocols.moves.AddPyMOLObserver(pose_pack_mover, True)
    pack_mover.apply(pose_pack_mover)

    print('Original Score: %s' % scorefxn(pose))
    print('Score after mutating 124, 123: %s' % scorefxn(pose_pack_mover))

    outFile.write('original energy: %.2f\n' % scorefxn(pose))
    outFile.write('energy after fast relax: %.2f\n' % scorefxn(pose_pack_mover))
    

def controlMut():
    from pyrosetta.toolbox import generate_resfile_from_pose
    generate_resfile_from_pose(pose, 'my.resfile')

def editResfile():
    import re
    global new_rfile
    fname = 'my'
    fext = 'resfile'
    rfile = open('.'.join([fname, fext]))
    text = rfile.read()
    rfile.close()
    res_to_mutate = [279, 275, 241, 332, 330, 359, 354, 276, 333, 254, 321, 277, 354]
    reg = r'( *(%s)\s*A\s*)(\w+)( *)' % ' | '.join(map(str,res_to_mutate))
    subStr = r'\1PIKAA AGFSTND\4'
    def subFun(m): return m.expand(subStr)         
    new_text = re.sub(reg,subFun,text)
    new_rfile = open('.'.join([fname + '_MUT2', fext]),'w+')
    new_rfile.write(new_text)
    new_rfile.close()


def packMutant():
    debugHeader('Packing Mutant')
    pose_pack_mover = Pose()
    pose_pack_mover.assign(pose)
    pose_pack_mover.pdb_info().name('defined_pack_mover')
  
    outFile.write('Beginning Packing of Mutant: {}\n'.format(now()))
    debugHeader("Writing File")
    mutResFileName = 'my_MUT2.resfile'
    outFile.write('Using Resfile: {}\n'.format(mutResFileName))
    debugHeader('File Being Written Again')

    debugHeader('Adding Observer')
    observer = rosetta.protocols.moves.AddPyMOLObserver(pose_pack_mover, True)

    kT = 1.0
    mc = MonteCarlo(pose_pack_mover, scorefxn, kT)
    
    for i in range(100):
        task = rosetta.core.pack.task.TaskFactory.create_packer_task(pose_pack_mover)
        rosetta.core.pack.task.parse_resfile(pose_pack_mover, task, mutResFileName)
        packer_task = rosetta.protocols.simple_moves.PackRotamersMover(scorefxn, task)

        min_mover = rosetta.protocols.simple_moves.MinMover()
        mm70150 = MoveMap()
        mm70150.set_bb_true_range(70, 150)
        min_mover.movemap(mm70150)
        min_mover.score_function(scorefxn)

        seq_mover = SequenceMover()
        seq_mover.add_mover(packer_task)
        seq_mover.add_mover(min_mover)

        trial_pack_min_mover = TrialMover(seq_mover, mc)

        trial_pack_min_mover.apply(pose_pack_mover)
        print('Score: %s' % scorefxn(pose_pack_mover))
        outFile.write('Iter: {0:d} --- Score: {1:.5f}\n'.format(i+1,scorefxn(pose_pack_mover)))
    
    outFile.write("Finished Packing of Mutant: {}\n".format(now()))
    # PyMOL> select my_ser, (resi 89+101+123+139+147) & defined_pack_mover

    
def main():
    # comment out sections you dont want to run...  for a given python
    # session, any function will run after all the previous functions
    # have run once

    debugHeader('Initialize')
    initialize()
    debugHeader('Load')
    load()
    debugHeader('Import PDB')
    importPdb()
    debugHeader('Send to PyMOL')
    sendtopymol()
    debugHeader('Show H Bonds')
    showhbonds()
    debugHeader('Find Minima')
    findmimima()

    # both of the fast relax sessions take a while...uncomment to run
    # debugHeader('Fast Relax 1')
    # fastrelax()
    debugHeader('Fast Relax 2')
    fastrelax2()

    debugHeader('Mutate Residues')
    mutataRes()
    debugHeader('COntrol Mutant')
    controlMut()
    debugHeader('Edit Resfile')
    editResfile()
    debugHeader('Packing Mutant')
    packMutant()
    debugHeader('Repeat Movers')
    repeatMovers()

main()
outFile.close()

