#!/usr/local/bin python3.6
# mainScript3_py3.py

import os
import re
import random
import pprint
from datetime import datetime

from pyrosetta import *
from pyrosetta.toolbox import generate_resfile_from_pose
from pyrosetta.rosetta.protocols.simple_moves import *
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.protocols.moves import AddPyMOLObserver
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task import parse_resfile


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

# FIXME - Integrate into ResfileBuilder Class
def makeSequenceGlobals():
    if madeGlobal('madeSeqGlobals'):
            return
    global pocketResNums, firstResNum, lastResNum,  oneIndexedRes,\
        allAa, aaGroupings, adnlAa, conservMutMap, liberalMutMap, madeSeqGlobals
    
    pocketResNums = [241, 247, 250, 251, 254, 255, 272, 273, 275, 276, 279, 280,
                     321, 325, 330, 332, 333, 334, 339, 343, 344, 354, 355, 358]
    firstResNum = 202
    lastResNum = 468
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
        for line in text.splitlines():
            f.writelines('[{}] >  {}\n'.format(now(2),line))

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


def createPyMolMover():
    pymol = PyMOLMover()
    return pymol


class CustomMover(pyrosetta.rosetta.protocols.moves.Mover):
    scorefxn = defaultScorefxn
    def __init__(self):
        pass
    
    def __str__(self):
        info = []
        for key,val in self.__dict__.items():
            infoStr = '{}: '.format(key)
            if isinstance(val,int):
                infoStr += '{:2d}'
            elif isinstance(val,float):
                infoStr += '{:7.1f}'
            else:
                infoStr += '{}'
            # print('Val: {} | InfoStr.format: {}'.format(val,infoStr.format(val)))
            infoStr = infoStr.format(val)
            info.append(infoStr)
        return ' | '.join(info)
    
    def get_name(self):
        return self.__class__.__name__


class FastRelaxMover(CustomMover):

    
    def __init__(self):
        pass

    def apply(self,pose):
        dprint('Beginning Fast Relax')
        fast_relax = FastRelax()
        fast_relax.set_scorefxn(self.scorefxn) 
        fast_relax.apply(pose)
    
    
class RepMinMover(CustomMover):


    def __init__(self):
        self.repeats = 10
        self.kT = 1.0
        self.bb_all = False
        self.bb_range = None
        self.bb_res = None
        self.chi_all = False
        self.chi_range = None
        self.chi_res = None
        
    def setMovemapBB(self,movemap):
        movemap.set_bb(False)
        if self.bb_res is not None:
            for i in self.bb_res:
                pose_i = pose.pdb_info().pdb2pose(i,'A')
                movemap.set_bb(pose_i,True)
        elif self.bb_range is not None:
            begin, end = self.bb_range
            movemap.set_bb_true_range(begin,end)
        elif self.bb_all:
            movemap.set_bb(True) 

    def setMovemapChi(self,movemap):
        movemap.set_chi(False)
        if self.chi_res is not None:
            for i in self.chi_res:
                pose_i = pose.pdb_info().pdb2pose(i,'A')
                movemap.set_chi(pose_i,True)
        elif self.chi_range is not None:
            begin, end = self.chi_range
            movemap.set_chi_true_range(begin,end)
        elif self.chi_all:
            movemap.set_chi(True) 
            
    def apply(self,pose):
        dprint('Perfoming A Minimization Loop')
        log( '`--> {}'.format(self))
        movemap = MoveMap()
        self.setMovemapBB(movemap)
        self.setMovemapChi(movemap)
        min_mv = MinMover()
        min_mv.movemap(movemap)
        min_mv.score_function(self.scorefxn)
        mc = MonteCarlo(pose,self.scorefxn,self.kT)
        trial_mv = TrialMover(min_mv,mc)
        rep_mv = RepeatMover(trial_mv,self.repeats)
        rep_mv.apply(pose)   

        
class SmallShearMover(CustomMover):


    def __init__(self):
        self.repeats = 50
        self.n_moves =5
        self.kT = 1.0
        self.angle = None
        self.bb_range = None
        
    def apply(self, pose):
        dprint('Running Small and Shear Movers'.format(repeats))
        log(' `--> {}'.format(self))

        min_mv = RepMinMover()
        min_mv.bb_range = self.bb_range
        min_mv.repeats = 2
        
        small_mv = SmallMover(movemap, self.kT, self.n_moves)
        shear_mv = ShearMover(movemap, self.kT, self.n_moves)
        if self.angle is not None:
            small_mv.angle_max('E', self.angle)
            small_mv.angle_max('L', self.angle)
            shear_mv.angle_max('E', self.angle)
            shear_mv.angle_max('L', self.angle)

        mc = MonteCarlo(pose, self.scorefxn, self.kT)

        seq_A = SequenceMover()
        seq_A.add_mover(small_mv)
        seq_A.add_mover(min_mv)
        trial_A = TrialMover(seq_A, mc)
        seq_B = SequenceMover()
        seq_B.add_mover(shear_mv)
        seq_B.add_mover(min_mv)
        trial_B = TrialMover(seq_B, mc)
        seq_total = SequenceMover()
        seq_total.add_mover(trial_A)
        seq_total.add_mover(trial_B)
        rep_mv = RepeatMover(seq_total, self.repeats)
        rep_mv.apply(pose)

        
class AnnealLoopMover(CustomMover):


    def __init__(self):
        self.cycles = 2
        self.kT = 10
        self.heat_time = 3
        self.anneal_time = 4
        self.angle_max = 4.0
        self.kT_max = 100.0
        self.seqRepeats = 5
        self.numSmallShearRepeats = 5
        self.kT_ratio = 20
        self.angle_ratio = 10
        self.bb_range = None
                
    def apply(self,pose):
        angles = ([self.angle_max] * self.heat_time
                  + [self.angle_max / self.angle_ratio] * self.anneal_time)
        kT2s = ([self.kT_max] * self.heat_time
                + [self.kT_max / self.kT_ratio] * self.anneal_time)
        cycleLen = len(angles)
        n1 = cycleLen* self.cycles
 
        ss_mv = SmallShearMover(self.seqRepeats,
                                self.numSmallShearRepeats,
                                kT,angle)
        if self.bb_range is not None:
            ss_mv.bb_range = self.bb_range
        mc = MonteCarlo(pose, scorefxn, self.kT)
        dprint('Beginning Anneal Loop')
        log(' `--> {}'.format(self))
        for i in range(n1):
            dprint('Beginning Loop # {:2d}/{:2d}'.format(i+1,n1))
            ind = i % cycleLen
            angle = angles[ind]
            kT2 = kT2s[ind]
            ss_mv.angle = angle
            ss_mv.kT = kT2
            ss_mv.apply(pose)
            mc.boltzmann(pose)
            printScore(pose,'Iteration # {:2d}'.format(i+1))


class MutantPackMover(CustomMover):


    def __init__(self):
        self.resfile = None
        pass

    def apply(self,pose):
        task = TaskFactory.create_packer_task(pose)
        parse_resfile(pose, task, self.resfile)
        pack_mv = PackRotamersMover(self.scorefxn, task)
        pack_mv.apply(pose)


class ResfileBuilder:

    
    resfile_dir = 'Resfiles'
    resfile_ext = 'resfile'
    rotamer_tag = 'NATAA'
    mutate_tag = 'PIKAA {}'
    chain = 'A'
    ligand_chain = 'X'
    line_fmt = '{:3d}  {}  {}\n'
    
    def __init__(self):
        self.filename = 'unnamed'
        self.packable_residues = []
        self.mutable_residues = []
        self.mut_liberal = False
        self.pose = []
        mkDir(self.resfile_dir)

    def getMutDict(self):
        if self.mut_liberal: return liberalMutMap
        else: return conservMutMap 

    def getFullFilename(self):
        return  self.filename + '.' + self.resfile_ext
        
    def getResfilePath(self):
        print('ResfileBuilder.getResfilePath: Joining {} and {}'.\
              format(self.resfile_dir, self.getFullFilename()))
        return os.path.join(self.resfile_dir,self.getFullFilename())

    def getFileHeader(self):
        text = ('# {}\n'.format(self.getFullFilename())
                + 'NATRO\nSTART\n\n'
                + '# Auto Generated Resfile by ResfileBuilder class\n' 
                + '# {}\n\n'.format(now()))
        return text
    
    def build(self):
        print('ResfileBuilder.build: filename --> {}'.format(self.filename))
        with open(self.getResfilePath(),'w') as rfile:
            rfile.write(self.getFileHeader())
            if len(self.packable_residues) > 0:
                rfile.write('# Packable Residues:\n')
            for res in self.packable_residues:
                line = self.line_fmt.format(res, self.chain, self.rotamer_tag)
                rfile.write(line)
            if len(self.mutable_residues) > 0:
                rfile.write('# Mutable Residues:\n')                
            for res in self.mutable_residues:
                line = self.line_fmt.format(res, self.chain, self.mutate_tag)
                poseNum = self.pose.pdb_info().pdb2pose(self.chain,res)
                aa = self.pose.residue(poseNum).name1()
                mutRes = mutDict[aa]
                line.format(''.join(mutRes))
                rfile.write(line)

    @classmethod
    def pocketRotamerResfile(cls):
        print('ResfileBuilder.pocketRotamerResfile: Creating Pocket Rotamer Resfile')
        name = 'pocket_rotamer'
        filepath = ResfileBuilder.resfilePath(name)
        if os.path.isfile(filepath):
            return filepath
        builder = cls()
        print('ResfileBuilder.pocketRotamerResfile: name --> {}'.format(name))
        builder.filename = name
        print('ResfileBuilder.pocketRotamerResfile: builder.filename --> {}'.\
              format(builder.filename))
        builder.packable_residues = pocketResNums
        builder.build()
        return builder.getResfilePath()

    @classmethod
    def fullRotamerResfile(cls):
        print('ResfileBuilder.fullRotamerResfile: Creating Full Rotamer Resfile')
        name = 'full_rotamer'
        filepath = ResfileBuilder.resfilePath(name)
        if os.path.isfile(filepath):
            return filepath
        builder = cls()
        builder.filename = name
        builder.packable_residues = range(firstResNum, lastResNum + 1)
        builder.build()
        return builder.getResfilePath()

    @classmethod
    def ligandRotamerResfile(cls):
        name = 'ligand_rotamer'
        filepath = ResfileBuilder.resfilePath(name)
        if os.path.isfile(filepath):
            return filepath
        builder = cls()
        builder.filename = name
        builder.packable_residues = [1]
        builder.chain = 'X'
        builder.build()
        return builder.getResfilePath()
        
    @classmethod
    def resfileFromDecoySpecs(cls,pose,residues,identifier,cycle)
        builder = cls()
        builder.mutable_residues = residues
        builder.filename = 'decoy{:02d}.{:02d}'.format(identifier,cycle)
        builder.pose = pose
        builder.build()
        return builder.getResfilePath()    
    
    @staticmethod
    def resfilePath(filename):
        suffix = '.' + ResfileBuilder.resfile_ext
        if not filename.endswith(suffix):
            filename += suffix
        return os.path.join(ResfileBuilder.resfile_dir,filename)
        

def rand_sample(n,lst):
    length = len(lst)
    temp_list = list(lst)
    sample = []
    for i in range(n):
        ind = random.randint(0, length - 1 - i)
        sample.append(temp_list[ind])
        del temp_list[ind]

    return sample,temp_list

class MutationLoopAssemblyMover(CustomMover):
    decoy_count = 0

    def __init__(self):
        self.mut_pattern = [[]]
        self.kT = 1.0
        self.annealLoop = AnnealLoopMover()
        self.identifier = MutationLoopMover.decoy_count
        self.decoy_count += 1
        
    def apply(self, pose=Pose()):
        for mut_residues in mut_pattern:
            for res in mut_residues:
                
          # for each decoy create a mover that:
          #   1. mutates the specified residues with resfile
          #   1. repack rotomers for all pocket residues, use movemap with MinMover

          
          #   1. do small/shear anneal loop
          self.annealLoop
          #   1. minimize pocket rotamers
          #   1. minimize backmobe
          # apply each specific mover to its decoy


        
        

def main():

    # Randomly sample 4 residues from the pocket residues
        # 6 times, so all that  pocket residues are used once
        # ^ do this n1 times
        pose = Pose()
        resPerDecoy = 4
        rand_samples = []
        n1 = 2
        for i in range(n1):
            remaining_res = list(pocketResNums)
            while len(remaining_res) >= 4:
                smp, remaining_res = rand_sample(resPerDecoy,remaining_res)
                rand_samples.append(smp)

        r_builder = ResfileBuilder()
        mut_files = []
        for sample in rand_sample:
            r_builder.mutable_residues = sample
            r_builder.filename = 'MUT_{}'.format('-'.join([str(i) for i in sample]))
            mut_files.append(r_builder.filename)
            r_builder.pose = pose
            r_builder.build()
    ndecoys = len(rand_samples)
    mkDir('Decoys')
    file_template = os.path.join('Decoys','output-{}')
    file_template.format(now(1))
    jd = PyJobDistributor(file_template,ndecoys,defaultScorefxn)
    while not jd.job_complete:
        new_pose.assign(pose)
        # Mover
        jd.output_decoy(pose)

if __name__ == '__main__':
    logBegin()
    main()
    logEnd()

 



