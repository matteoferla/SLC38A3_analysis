# ED and Membrane

With the hindsight knowledge that the extracellular fluff was causing low similarity, 
threading EM-fluff–stripped SLC38A3 against electron-density in membrane minimised PDB:6C08 
with `pyrosetta.rosetta.protocols.comparative_modeling.ThreadingMover` might have been better for the minimation step.
This is because membrane protein can thread weirdly in I-Tasser and Phyre.

## Commands

The following command add a membrane vector for a OPM-aligned structure.
Pyrosetta only cares 
OMP files have `REMARK      1/2 of bilayer thickness:   14.3` and a `CRYST` entry
These are not important.
Nor are the DUM atoms, elements O/N ± 14 Å on z making two leaves in a 2x2 Å grid.
What is important is that the z=0 plane is the middle of the membrane and
the z=±14.5 Å planes are the upper/lower boundaries.

    addmem = pyrosetta.rosetta.protocols.membrane.AddMembraneMover("from_structure")
    addmem.apply(pose)
    
The alternative, a spanfile, works only if the structure is OPM-aligned, even if it is the spanfile generated from the above command (as `out.span`)
Otherwise it segfaults

    addmem = pyrosetta.rosetta.protocols.membrane.AddMembraneMover('out.span')
    
This is a channel, so adding a pore is beneficial

    pyrosetta.rosetta.protocols.membrane.AqueousPoreFinder().apply(pose)

The trick is then to align it within Pyrosetta
    
    ali = pyrosetta.rosetta.protocols.simple_moves.AlignChainMover()
    ali.pose(map_pose)
    ali.source_chain(1)  #whole structure would give ERROR: Assertion `count == ref_coords.size()` failed.
    ali.target_chain(1)
    ali.apply(mem_pose)
    
In the presence of ligands, loading the polymer, adding the membrane and then adding the ligands allows it work.

    with pymol2.PyMOL() as pymol:
        pymol.cmd.load(filename, 'S38A3')
        map_pdbblock = pymol.cmd.get_pdbstr()
        pymol.cmd.load(opm_filename, 'OPM')
        pymol.cmd.align('S38A3', 'OPM')
        pymol.cmd.remove('resn MEM') # make sure MEM is not there!
        pymol.cmd.delete('OPM')
        holo_pdbblock = pymol.cmd.get_pdbstr('polymer')
        lig_pdbblock = pymol.cmd.get_pdbstr('not polymer')
    pose = pyrosetta.Pose()
    if params_paths:
        pyrosetta.generate_nonstandard_residue_set(pose, params_paths)
    pyrosetta.rosetta.core.import_pose.pose_from_pdbstring(pose, holo_pdbblock)

    addmem = pyrosetta.rosetta.protocols.membrane.AddMembraneMover('from_structure')
    addmem.apply(pose)
    pyrosetta.rosetta.protocols.membrane.AqueousPoreFinder().apply(pose)
    ligpose = pyrosetta.Pose()
    pyrosetta.generate_nonstandard_residue_set(ligpose, params_paths)
    pyrosetta.rosetta.core.import_pose.pose_from_pdbstring(ligpose, lig_pdbblock)
    pose.append_pose_by_jump(ligpose, pose.total_residue())
    return pose

    
The [electron density minimisation described elsewhere](http://blog.matteoferla.com/2020/04/how-to-set-up-electron-density.html) 
works fine with the membrane.
The next issue is to prevent the protein from boobing out of the membrane, when not constrained to an electron density.
As a consequence, a series of `CoordinateConstraint` restraints to the centre atom of the MEM vector to the residues (pose numbered)
2, 203, 318, 349, 383, 389 were added, first with &times;100, then &times;1 weight.

NB. the `pdb_info` gets blanked on `.dump_pdb(...)`.
