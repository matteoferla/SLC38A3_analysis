import pymol2, pyrosetta


def load_mem_pose(filename, opm_filename, params_paths=None):
    with pymol2.PyMOL() as pymol:
        #     pymol.cmd.load('SLC38A3_phyre.pdb', 'S38A3')
        pymol.cmd.load(filename, 'S38A3')
        # pymol.cmd.remove('resi 1-68 or resi 245-282 or resi 494-504')
        #     pymol.cmd.remove('resi 138-207 or resi 318-392')
        map_pdbblock = pymol.cmd.get_pdbstr()
        pymol.cmd.load(opm_filename, 'OPM')
        pymol.cmd.align('S38A3', 'OPM')
        pymol.cmd.remove('resn MEM')  # make sure its not there!
        pymol.cmd.delete('OPM')
        holo_pdbblock = pymol.cmd.get_pdbstr('polymer')
        lig_pdbblock = pymol.cmd.get_pdbstr('not polymer')
    #     pymol.cmd.save('test.pdb')
    pose = pyrosetta.Pose()
    if params_paths:
        pyrosetta.generate_nonstandard_residue_set(pose, params_paths)
    # pyrosetta.rosetta.core.import_pose.pose_from_file(pose, 'phyre.gln_Na.PO4.loop_r.pdb')
    pyrosetta.rosetta.core.import_pose.pose_from_pdbstring(pose, holo_pdbblock)

    addmem = pyrosetta.rosetta.protocols.membrane.AddMembraneMover('from_structure')
    addmem.apply(pose)
    pyrosetta.rosetta.protocols.membrane.AqueousPoreFinder().apply(pose)
    ligpose = pyrosetta.Pose()
    pyrosetta.generate_nonstandard_residue_set(ligpose, params_paths)
    # pyrosetta.rosetta.core.import_pose.pose_from_file(pose, 'phyre.gln_Na.PO4.loop_r.pdb')
    pyrosetta.rosetta.core.import_pose.pose_from_pdbstring(ligpose, lig_pdbblock)
    pose.append_pose_by_jump(ligpose, pose.total_residue())
    return pose
