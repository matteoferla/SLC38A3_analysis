# ED without membrane
# ======== init pyrosetta ==============================================================================================
import pyrosetta
from .init_boilerplate import make_option_string
pyrosetta.init(extra_options=make_option_string(no_optH=False,
                                                ex1=None,
                                                ex2=None,
                                                mute='all',
                                                ignore_unrecognized_res=True,
                                                load_PDB_components=False,
                                                ignore_waters=False)
               )

# ======== Pose ========================================================================================================
filename = 'SLC38A3_phyre.pdb'
import pymol2

with pymol2.PyMOL() as pymol:
    pymol.cmd.load(filename, 'S38A3')
    pymol.cmd.remove('resi 1-68 or resi 245-282 or resi 494-504')
    pymol.cmd.alter('*', 'chain="A"')
    pymol.cmd.sort()
    pymol.cmd.fetch('6C08')
    pymol.cmd.align('S38A3', '6C08 and chain C')
    pymol.cmd.delete('6C08')
    pdbblock = pymol.cmd.get_pdbstr()

pose = pyrosetta.Pose()
pyrosetta.rosetta.core.import_pose.pose_from_pdbstring(pose, pdbblock)

# ======== scorefxn ====================================================================================================

scorefxnED = pyrosetta.get_fa_scorefxn()  # ref2015 --> franklin2019
ED = pyrosetta.rosetta.core.scoring.electron_density.getDensityMap('6c08.ccp4')  # from PDBe
sdsm = pyrosetta.rosetta.protocols.electron_density.SetupForDensityScoringMover()
sdsm.apply(pose)
elec_dens_fast = pyrosetta.rosetta.core.scoring.ScoreType.elec_dens_fast
scorefxnED.set_weight(elec_dens_fast, 30)

# ======== relax =======================================================================================================

for w in (30, 20, 10):
    scorefxnED.set_weight(elec_dens_fast, w)
    relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxnED, 5)
    relax.apply(pose)
    print('Match to ED', ED.matchPose(pose))
    print('ED weight', scorefxnED.weights()[elec_dens_fast])
    print('score', scorefxnED(pose))
reg_scorefxn = pyrosetta.get_fa_scorefxn()
print('How good is the density match:', ED.matchPose(pose))
print(reg_scorefxn.get_name(), reg_scorefxn(pose))

