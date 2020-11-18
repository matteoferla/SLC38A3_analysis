import pymol2

def fix_pdb(in_filename, out_filename=None):
    with pymol2.PyMOL() as pymol:
        pymol.cmd.load(in_filename, 'S38A3')
        pymol.cmd.alter('chain A', 'resv+=69-1')
        pymol.cmd.sort()
        pymol.cmd.alter('chain A and resi 245-9999', 'resv+=283-245')
        pymol.cmd.sort()
        if out_filename is not None:
            pymol.cmd.save(out_filename)
        return pymol.cmd.get_pdbstr()
