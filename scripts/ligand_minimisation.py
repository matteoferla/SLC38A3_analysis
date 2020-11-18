"""
There is a slight problem with phosphates on the membrane boundary, namely the part that is in the implicit surface wants to jump out.
Consequently for a final energy minimisation several loose constraints are needed to prevent them from being cannonballed out of the membrane.

"""

# ED without membrane
# ======== init pyrosetta ==============================================================================================
import pyrosetta
from .init_boilerplate import make_option_string
from .load_mem_pose import load_mem_pose

pyrosetta.init(extra_options=make_option_string(no_optH=False,
                                                ex1=None,
                                                ex2=None,
                                                mute='all',
                                                ignore_unrecognized_res=True,
                                                load_PDB_components=False,
                                                ignore_waters=False)
               )

# ======== Pose ========================================================================================================

filename = 'phyre.gln_Na.PO4.numb.pdb'
opm_filename = 'OPM_6c08.pdb'
params_paths = pyrosetta.rosetta.utility.vector1_string(1)
params_paths[1] = 'PO4.params'
pose = load_mem_pose(filename, opm_filename, params_paths)
scorefxn = pyrosetta.create_score_function('franklin2019')

# ======== Common ======================================================================================================

# constraints require atomID, which are a bit clunky to select so here's a function to select them.
def get_AtomID(chain:str, resi:int, atomname:str) -> pyrosetta.rosetta.core.id.AtomID:
    r = pose.pdb_info().pdb2pose(res=resi, chain=chain)
    assert r != 0, f'{resi}:{chain} is absent'
    residue = pose.residue(r)
    return pyrosetta.rosetta.core.id.AtomID(atomno_in=residue.atom_index(atomname), rsd_in=r)


HarmonicFunc = pyrosetta.rosetta.core.scoring.func.HarmonicFunc
# see pyrosetta.rosetta.core.scoring.func for the others!

# first lets make a list...
cons = []

# ======== coordinate constaint ========================================================================================
# the protein may dive out of the membrane

constrained = [  [203, -22.254404067993164, -7.722462177276611, -2.72951340675354],
                 [318, 7.041024208068848, 11.991669654846191, -18.32436752319336],
                 [349, 4.768013000488281, 5.710573196411133, 18.04562759399414],
                 [383, -12.9513521194458, -5.184802055358887, -10.952046394348145],
                 [389, -9.103133201599121, -14.593981742858887, -11.054950714111328]
              ]
nsv = pyrosetta.rosetta.core.select.residue_selector.ResidueNameSelector('MEM').apply(pose)
r = pyrosetta.rosetta.core.select.residue_selector.ResidueVector(nsv)[1]
CNTR = pyrosetta.rosetta.core.id.AtomID(atomno_in=pose.residue(r).atom_index('CNTR'),
                                        rsd_in=r)

# Specify the atoms you want...
CoordinateConstraint = pyrosetta.rosetta.core.scoring.constraints.CoordinateConstraint
for resi, x, y, z in constrained:
    cons.append(CoordinateConstraint(get_AtomID('A', resi, 'CA'),
                                     CNTR,
                                     pyrosetta.rosetta.numeric.xyzVector_double_t(x, y, z),
                                     HarmonicFunc(x0_in=0, sd_in=0.5)
                                    )
               )

# ======== residue constaint ===========================================================================================

## GLN as per PDB

cons.append( pyrosetta.rosetta.core.scoring.constraints.AtomPairConstraint(get_AtomID('A', 160, 'OH'),
                                                                    get_AtomID('B', 388, 'O'),
                                                                    HarmonicFunc(x0_in=3.4, sd_in=0.2))

           )

cons.append( pyrosetta.rosetta.core.scoring.constraints.AtomPairConstraint(get_AtomID('A', 80, 'O'),
                                                                    get_AtomID('B', 388,'NE2'),
                                                                    HarmonicFunc(x0_in=3.4, sd_in=0.2))

           )

cons.append( pyrosetta.rosetta.core.scoring.constraints.AtomPairConstraint(get_AtomID('A', 82, 'OG'),
                                                                    get_AtomID('B', 388, 'N'),
                                                                    HarmonicFunc(x0_in=3.4, sd_in=0.2))

           )

## Phospho
phosphoarg_pdbs = {391: [208,390,389], 393: [135, 454]}

for phospho in phosphoarg_pdbs:
    for arg in phosphoarg_pdbs[phospho]:
        cons.append( pyrosetta.rosetta.core.scoring.constraints.AtomPairConstraint(get_AtomID('A', arg, 'CZ'),
                                                                        get_AtomID('E', phospho, 'P'),
                                                                        HarmonicFunc(x0_in=4.3, sd_in=0.4))

               )
cons.append( pyrosetta.rosetta.core.scoring.constraints.AtomPairConstraint(get_AtomID('A', 389, '1HH2'),
                                                                get_AtomID('E', 391, 'O1'),
                                                                HarmonicFunc(x0_in=2.3, sd_in=0.2))

       )

## For the sake of principle... hydrogen of PO4 membrane-wards.
cons.append( pyrosetta.rosetta.core.scoring.constraints.AtomPairConstraint(get_AtomID('A', 142, 'CE'),
                                                                get_AtomID('E', 393, 'P'),
                                                                HarmonicFunc(x0_in=4.3, sd_in=1))

       )

## Prevent wandering C-terminus


#      You clicked /combo//A/VAL`198/CA -> (pk1)
#  You clicked /combo//A/ASN`419/CA -> (pk2)
cons.append( pyrosetta.rosetta.core.scoring.constraints.AtomPairConstraint(get_AtomID('A', 198, 'CA'),
                                                                get_AtomID('A', 419, 'CA'),
                                                                HarmonicFunc(x0_in=6.798631191253662, sd_in=0.01))

       )

#'/model1//A/TRP`493/CA', '/model1//A/ASN`283/CA')
#
cons.append( pyrosetta.rosetta.core.scoring.constraints.AtomPairConstraint(get_AtomID('A', 493, 'CA'),
                                                                get_AtomID('A', 283, 'CA'),
                                                                HarmonicFunc(x0_in=18.57731819152832, sd_in=0.01))

       )

# PyMOL>print cmd.distance('/model1//A/PRO`244/CA', '/model1//A/ASN`283/CA')
#
cons.append( pyrosetta.rosetta.core.scoring.constraints.AtomPairConstraint(get_AtomID('A', 244, 'CA'),
                                                                get_AtomID('A', 283, 'CA'),
                                                                HarmonicFunc(x0_in=18.179616928100586, sd_in=0.01))

       )

# ======== Common ======================================================================================================

# lets make the list into a vector
cl = pyrosetta.rosetta.utility.vector1_std_shared_ptr_const_core_scoring_constraints_Constraint_t()
cl.extend(cons)

# lets make the vector into a constraint set
cs = pyrosetta.rosetta.core.scoring.constraints.ConstraintSet()
cs.add_constraints(cl)

setup = pyrosetta.rosetta.protocols.constraint_movers.ConstraintSetMover()
setup.constraint_set(cs)
setup.apply(pose)

stm = pyrosetta.rosetta.core.scoring.ScoreTypeManager()
scorefxn.set_weight(stm.score_type_from_name("atom_pair_constraint"), 2)
scorefxn.set_weight(stm.score_type_from_name("coordinate_constraint"), 5)

# ======== Relax =======================================================================================================

relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, 5)
movemap = pyrosetta.MoveMap()
mem_sele = pyrosetta.rosetta.core.select.residue_selector.ResidueNameSelector('MEM')
not_mem_sele = pyrosetta.rosetta.core.select.residue_selector.NotResidueSelector(mem_sele)
nmv = not_mem_sele.apply(pose)
movemap.set_bb(allow_bb=nmv)
movemap.set_chi(allow_chi=nmv)
movemap.set_jump(True)
# relax.set_native_pose(pose.clone())
# relax.constrain_relax_to_native_coords(True)
relax.set_movemap(movemap)
# relax.constrain_relax_to_start_coords(True)


pose.dump_pdb('test-pre.pdb')

relax.apply(pose)

pose.dump_pdb('test-post.pdb')


