# This has several modification relative to previous code. Namely, `Mutation` has been split out.

import pyrosetta
import re, os, csv
from typing import Optional, List, Dict, Union


class Mutation:
    """
    A mutation is an object that has all the details of the mutation.
    A variant, as interpreted in ``Model.score_mutations`` is a pose with a mutation.
    """
    _name3 = {'A': 'ALA',
              'C': 'CYS',
              'D': 'ASP',
              'E': 'GLU',
              'F': 'PHE',
              'G': 'GLY',
              'H': 'HIS',
              'I': 'ILE',
              'L': 'LEU',
              'K': 'LYS',
              'M': 'MET',
              'N': 'ASN',
              'P': 'PRO',
              'Q': 'GLN',
              'R': 'ARG',
              'S': 'SER',
              'T': 'THR',
              'V': 'VAL',
              'W': 'TRP',
              'Y': 'TYR'}

    def __init__(self, mutation: str, chain: str, pose: pyrosetta.Pose):
        self.mutation = self.parse_mutation(mutation)
        rex = re.match('(\w)(\d+)(\w)', self.mutation)
        self.pdb_resi = int(rex.group(2))
        self.chain = chain
        self.from_resn1 = rex.group(1)
        self.from_resn3 = self._name3[rex.group(1)]
        self.to_resn1 = rex.group(3)
        self.to_resn3 = self._name3[rex.group(3)]
        pose2pdb = pose.pdb_info().pdb2pose
        self.pose_resi = pose2pdb(res=self.pdb_resi, chain=self.chain)
        if self.pose_resi != 0:
            self.pose_residue = pose.residue(self.pose_resi)
            self.pose_resn1 = self.pose_residue.name1()
            self.pose_resn3 = self.pose_residue.name3()
        else:
            self.pose_residue = None
            self.pose_resn1 = None
            self.pose_resn3 = None

    def parse_mutation(self, mutation):
        if mutation[:2] == 'p.':
            mutation = mutation.replace('p.', '')
        if mutation[1].isdigit():
            return mutation
        else:
            value2key = lambda value: list(self._name3.keys())[list(self._name3.values()).index(value.upper())]
            return value2key(mutation[:3]) + mutation[3:-3] + value2key(mutation[-3:])

    def is_valid(self):
        return self.pose_resn1 == self.from_resn1

    def assert_valid(self):
        assert self.is_valid(), f'residue {self.pose_resi}(pose)/{self.pdb_resi}(pdb) ' + \
                                f'is a {self.pose_resn3}, not a {self.from_resn3}'

    def __str__(self):
        return self.mutation


class Variant:
    """
    Copy pasted from PI4KA <- GNB2 <- SnoopCatcher
    """
    strict_about_starting_residue = True
    scorefxn = pyrosetta.get_fa_scorefxn()

    def __init__(self, pose: pyrosetta.Pose):
        self.pose = pose

    @classmethod
    def from_file(cls, filename: str, params_filenames: Optional[List[str]] = None):
        return cls(pose=cls.load_pose_from_file(filename, params_filenames))

    @classmethod
    def load_pose_from_file(self, filename: str, params_filenames: Optional[List[str]] = None) -> pyrosetta.Pose:
        """
        Loads a pose from filename with the params in the params_folder
        :param filename:
        :param params_filenames:
        :return:
        """
        pose = pyrosetta.Pose()
        if params_filenames:
            params_paths = pyrosetta.rosetta.utility.vector1_string()
            params_paths.extend(params_filenames)
            pyrosetta.generate_nonstandard_residue_set(pose, params_paths)
        pyrosetta.rosetta.core.import_pose.pose_from_file(pose, filename)
        return pose

    def relax_around_mover(self,
                           pose: pyrosetta.Pose,
                           mutation=None,
                           resi: int = None, chain: str = None, cycles=5, distance=5,
                           cartesian=False, own_chain_only=False) -> None:
        """
        Relaxes pose ``distance`` around resi:chain or mutation

        :param resi: PDB residue number.
        :param chain:
        :param pose:
        :param cycles: of relax (3 quick, 15 thorough)
        :param distance:
        :param cartesian:
        :return:
        """
        if mutation is None and resi is None:
            raise ValueError('mutation or resi+chain required')
        elif mutation is not None:
            resi = mutation.pose_resi
            chain = None
        else:
            pass
        if pose is None:
            pose = self.pose
        movemap = pyrosetta.MoveMap()
        ####
        n = self.get_neighbour_vector(pose=pose, resi=resi, chain=chain, distance=distance,
                                      own_chain_only=own_chain_only)
        # print(pyrosetta.rosetta.core.select.residue_selector.ResidueVector(n))
        movemap.set_bb(False)
        movemap.set_bb(allow_bb=n)
        movemap.set_chi(False)
        movemap.set_chi(allow_chi=n)
        movemap.set_jump(False)
        relax = pyrosetta.rosetta.protocols.relax.FastRelax(self.scorefxn, cycles)
        relax.set_movemap(movemap)
        relax.set_movemap_disables_packing_of_fixed_chi_positions(True)
        relax.cartesian(cartesian)
        relax.apply(pose)

    def get_neighbour_vector(self, pose: pyrosetta.Pose, resi: int, chain: str, distance: int,
                             include_focus_in_subset: bool = True,
                             own_chain_only: bool = False) -> pyrosetta.rosetta.utility.vector1_bool:
        resi_sele = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
        if chain is None:  # pose numbering.
            resi_sele.set_index(resi)
        else:
            resi_sele.set_index(pose.pdb_info().pdb2pose(chain=chain, res=resi))
        NeighborhoodResidueSelector = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector
        neigh_sele = NeighborhoodResidueSelector(resi_sele, distance=distance,
                                                 include_focus_in_subset=include_focus_in_subset)
        if own_chain_only and chain is not None:
            chain_sele = pyrosetta.rosetta.core.select.residue_selector.ChainSelector(chain)
            and_sele = pyrosetta.rosetta.core.select.residue_selector.AndResidueSelector(neigh_sele, chain_sele)
            return and_sele.apply(pose)
        else:
            return neigh_sele.apply(pose)

    def parse_mutation(self, mutation: Union[str, Mutation], chain, pose: pyrosetta.Pose = None):
        if pose is None:
            pose = self.pose
        if isinstance(mutation, str):
            mutant = Mutation(mutation, chain, pose)
        elif isinstance(mutation, Mutation):
            mutant = mutation
        else:
            raise TypeError(f'Does not accept mutation of type {mutation.__class__.__name__}')
        if mutant.pose_resi == 0:
            raise ValueError('Not in pose')
        if self.strict_about_starting_residue:
            mutant.assert_valid()
        return mutant

    def make_mutant(self,
                    pose: pyrosetta.Pose,
                    mutation: Union[str, Mutation],
                    chain='A',
                    distance: int = 10,
                    cycles: int = 5
                    ) -> pyrosetta.Pose:
        """
        Make a point mutant (``A23D``).
        :param pose: pose
        :param mutation:
        :param chain:
        :return:
        """
        if pose is None:
            mutant = self.pose.clone()
        else:
            mutant = pose.clone()
        if isinstance(mutation, str):
            mutation = Mutation(mutation, chain, mutant)
        MutateResidue = pyrosetta.rosetta.protocols.simple_moves.MutateResidue
        MutateResidue(target=mutation.pose_resi, new_res=mutation.to_resn3).apply(mutant)
        self.relax_around_mover(mutant,
                                mutation=mutation,
                                distance=distance,
                                cycles=cycles,
                                own_chain_only=False)
        return mutant

    def score_interface(self, pose: pyrosetta.Pose, interface: str) -> Dict[str, float]:
        if pose is None:
            pose = self.pose
        assert self.has_interface(pose, interface), f'There is no {interface}'
        ia = pyrosetta.rosetta.protocols.analysis.InterfaceAnalyzerMover(interface)
        ia.apply(pose)
        return {'complex_energy': ia.get_complex_energy(),
                'separated_interface_energy': ia.get_separated_interface_energy(),
                'complexed_sasa': ia.get_complexed_sasa(),
                'crossterm_interface_energy': ia.get_crossterm_interface_energy(),
                'interface_dG': ia.get_interface_dG(),
                'interface_delta_sasa': ia.get_interface_delta_sasa()}

    def has_interface(self, pose: pyrosetta.Pose, interface: str) -> bool:
        if pose is None:
            pose = self.pose
        pose2pdb = pose.pdb_info().pose2pdb
        have_chains = {pose2pdb(r).split()[1] for r in range(1, pose.total_residue() + 1)}
        want_chains = set(interface.replace('_', ''))
        return have_chains == want_chains

    def has_residue(self, pose: pyrosetta.Pose, resi: int, chain: str) -> bool:
        if pose is None:
            pose = self.pose
        pdb2pose = pose.pdb_info().pdb2pose
        r = pdb2pose(res=resi, chain=chain)
        return r != 0

    def vector2list(self, vector: pyrosetta.rosetta.utility.vector1_bool) -> pyrosetta.rosetta.std.list_unsigned_long_t:
        rv = pyrosetta.rosetta.core.select.residue_selector.ResidueVector(vector)
        x = pyrosetta.rosetta.std.list_unsigned_long_t()
        assert len(rv) > 0, 'Vector is empty!'
        for w in rv:
            x.append(w)
        return x

    def CA_RMSD(self, poseA: pyrosetta.Pose, poseB: pyrosetta.Pose, resi: int, chain: str, distance: int) -> float:
        n = self.get_neighbour_vector(pose=poseA, resi=resi, chain=chain, distance=distance)
        residues = self.vector2list(n)
        return pyrosetta.rosetta.core.scoring.CA_rmsd(poseA, poseB, residues)

    def FA_RMSD(self, poseA: pyrosetta.Pose, poseB: pyrosetta.Pose, resi: int, chain: str, distance: int) -> float:
        n = self.get_neighbour_vector(pose=poseA, resi=resi, chain=chain, distance=distance,
                                      include_focus_in_subset=False)
        residues = self.vector2list(n)
        # pyrosetta.rosetta.core.scoring.automorphic_rmsd(residueA, residueB, False)
        return pyrosetta.rosetta.core.scoring.all_atom_rmsd(poseA, poseB, residues)

    def does_contain(self, mutation: Union[Mutation, str], chain: Optional[str] = None) -> bool:
        assert isinstance(mutation, Mutation) or chain is not None, 'mutation as str requires chain.'
        if isinstance(mutation, str):
            mutation = Mutation(mutation, chain, self.pose)
        if mutation.pose_resi == 0:
            return False
        if mutation.pose_resn1 not in mutation._name3.keys():
            return False
        else:
            return True

    def score_mutations(self,
                        mutations,
                        chain='A',
                        interfaces=(),
                        modelname='holo',
                        preminimise=False,
                        verbose=True,
                        distance=10,
                        cycles=5):
        if not os.path.exists('variants'):
            os.mkdir('variants')
        with open(modelname + '_scores.csv', 'w') as w:
            fieldnames = ['model',
                          'mutation',
                          'complex_ddG',
                          'complex_native_dG',
                          'complex_mutant_dG',
                          'FA_RMSD',
                          'CA_RMSD'] + [f'{interface_name}_{suffix}' for suffix in ('interface_native_dG',
                                                                                    'interface_mutant_dG',
                                                                                    'interface_ddG')
                                        for interface_name, interface_scheme in interfaces]

            out = csv.DictWriter(w, fieldnames=fieldnames)
            out.writeheader()
            ## wt
            ref_interface_dG = {}
            if not preminimise:
                n = self.scorefxn(self.pose)
                for interface_name, interface_scheme in interfaces:
                    ref_interface_dG[interface_name] = self.score_interface(self.pose, interface_scheme)['interface_dG']
            else:
                pass  # calculated for each.
            ## muts
            for mut in mutations:
                try:
                    mutation = self.parse_mutation(mut, chain)
                    if verbose:
                        print(mutation)
                    if not self.does_contain(mutation):
                        if verbose:
                            print('Absent')
                        continue
                    if preminimise:
                        premutant = self.pose.clone()
                        self.relax_around_mover(premutant,
                                                mutation=mutation,
                                                distance=distance,
                                                cycles=cycles)
                        n = self.scorefxn(premutant)
                    else:
                        premutant = self.pose
                    variant = self.make_mutant(premutant,
                                               mutation=mutation,
                                               distance=distance,
                                               cycles=cycles)
                    variant.dump_pdb(f'variants/{modelname}.{mutation}.pdb')
                    m = self.scorefxn(variant)
                    data = {'model': modelname,
                            'mutation': str(mutation),
                            'complex_ddG': m - n,
                            'complex_native_dG': n,
                            'complex_mutant_dG': m,
                            'FA_RMSD': self.FA_RMSD(self.pose,
                                                    variant,
                                                    resi=mutation.pose_resi,
                                                    chain=None,
                                                    distance=distance),
                            'CA_RMSD': self.CA_RMSD(self.pose,
                                                    variant,
                                                    resi=mutation.pose_resi,
                                                    chain=None,
                                                    distance=distance)
                            }
                    for interface_name, interface_scheme in interfaces:
                        if self.has_interface(variant, interface_scheme):
                            if preminimise:
                                ref_interface_dG[interface_name] = self.score_interface(premutant, interface_scheme)[
                                    'interface_dG']
                            print(f'{interface_name} ({interface_scheme}) applicable to {modelname}')
                            i = self.score_interface(variant, interface_scheme)['interface_dG']
                        else:
                            print(f'{interface_name} ({interface_scheme}) not applicable to {modelname}')
                            i = float('nan')
                        data[f'{interface_name}_interface_native_dG'] = ref_interface_dG[interface_name]
                        data[f'{interface_name}_interface_mutant_dG'] = i
                        data[f'{interface_name}_interface_ddG'] = i - ref_interface_dG[interface_name]
                except Exception as error:
                    msg = f"{error.__class__.__name__}: {error}"
                    print(msg)
                    data = {'model': modelname,
                            'mutation': str(mut),
                            'complex_ddG': msg
                            }
                out.writerow(data)
