import sys
import pandas as pd
import numpy as np
from .io import import_para
from .io import import_input
from .io import import_config
from .io import import_dumppos
from .io import import_dumpbond
from .io import import_xyz
from .io import export_input
from .io import export_dumppos
from .io import export_xyz
from .analysis import bond_analysis
from .analysis import molecule_analysis
from .analysis import atom_analysis


class SimulationDat():
    def __init__(self):
        self.atoms = None
        self.cell = [None] * 3
        self.bondorder_list = None
        self.bondorder_connect_list = None
        self.connect_list = None
        self.connect_list_cut_off = None

        # atom set which exist in the system
        self.atom_type_set = set()

        # variables for para.rd
        self.atom_symbol_to_type = None
        self.atom_type_to_symbol = None
        self.atom_type_to_mass = None

        # variables for config.rd
        self.mpigrid = [1] * 3
        self.ompgrid = [1] * 3
        self.cutoff = 2.5
        self.margin = 0.3
        self.total_step = 0
        self.file_step = 1000
        self.flagdecomp = False
        self.flagconnect = False

        # calclation infomation in input files
        self.fix_info = []
        self.move_info = []
        self.press_info = []
        self.sumforce_info = []
        self.thermofree_info = []
        self.wall_info = []

    # IMPORT
    def import_input(self, ifn: str) -> None:
        importer = import_input.ImportInput()
        importer.import_file(self, ifn)

    def import_config(self, ifn: str) -> None:
        importer = import_config.ImportConfig()
        importer.import_file(self, ifn)

    def import_dumppos(self, ifn: str) -> None:
        importer = import_dumppos.ImportDumppos()
        importer.import_file(self, ifn)

    def import_dumpbond(self, ifn: str) -> None:
        importer = import_dumpbond.ImportDumpbond()
        importer.import_file(self, ifn)

    def import_xyz(self, ifn: str) -> None:
        importer = import_xyz.ImportXyz()
        importer.import_file(self, ifn)

    def import_para(self, ifn: str) -> None:
        importer = import_para.ImportPara()
        importer.import_file(self, ifn)

    # EXPORT
    def export_input(self, ofn: str) -> None:
        importer = to_input.ExportInput()
        importer.to_file(self, ofn)

    def export_dumppos(self, ofn: str, time_step=None, out_columns=None) -> None:
        importer = to_dumppos.ExportDumppos()
        importer.to_file(self, ofn, time_step, out_columns)

    def export_xyz(self, ofn: str, out_columns=None, structure_name=None) -> None:
        importer = to_xyz.ExportXyz()
        importer.to_file(self, ofn, out_columns, structure_name)

    # METHODS\
    def get_total_atoms(self):
        if self.atoms is not None:
            return len(self.atoms)
        elif self.bondorder_list is not None:
            return len(self.bondorder_list)
        else:
            print('Import file first')
            sys.exit(-1)

    def set_decomp(self):
        clx = float(1.0 / (self.ompgrid[0] * self.mpigrid[0]))
        cly = float(1.0 / (self.ompgrid[1] * self.mpigrid[1]))
        clz = float(1.0 / (self.ompgrid[2] * self.mpigrid[2]))
        self.decompx = [clx] * self.ompgrid[0] * self.mpigrid[0]
        self.decompy = [cly] * self.ompgrid[1] * self.mpigrid[1]
        self.decompz = [clz] * self.ompgrid[2] * self.mpigrid[2]

    def wrap_particles(self):
        self.atoms[['x', 'y', 'z']] %= self.cell

    def get_connect_list(self, cut_off):
        if self.connect_list_cut_off != cut_off or self.connect_list is None:
            self.__create_connect_list(cut_off)
        self.connect_list_cut_off = cut_off
        return self.connect_list

    def __create_connect_list(self, cut_off):
        if cut_off is None:
            print('cut_off is not defined')
            sys.exit(-1)
        if self.bondorder_list is None:
            print('bondorder_list is not defined')
            print('Import dumppos first')
            sys.exit(-1)
        self.connect_list = [[] for _ in range(self.get_total_atoms())]
        for atom_idx, (neibour_idxs, bondorder_list) in enumerate(zip(self.bondorder_connect_list, self.bondorder_list)):
            for neibour_idx, bond_l in zip(neibour_idxs, bondorder_list):
                if bond_l[-1] >= cut_off:
                    self.connect_list[atom_idx].append(neibour_idx)

    def replicate_atoms(self, replicate_directions=[1, 1, 1]):
        shift = self.cell
        shifted_atoms_list = [self.atoms]

        for x_idx in range(replicate_directions[0]):
            for y_idx in range(replicate_directions[1]):
                for z_idx in range(replicate_directions[2]):
                    if x_idx == 0 and y_idx == 0 and z_idx == 0:
                        continue
                    shifted_atoms = self.atoms.copy()
                    shifted_atoms[['x', 'y', 'z']] += shift * \
                        np.array([x_idx, y_idx, z_idx])
                    shifted_atoms_list.append(shifted_atoms)
        self.atoms = pd.concat(shifted_atoms_list)
        self.atoms.reset_index(drop=True, inplace=True)
        for dim in range(3):
            self.cell[dim] *= replicate_directions[dim]

    def concat_atoms(self, outer_sdat):
        self.atoms = pd.concat([self.atoms, outer_sdat.atoms])
        self.atoms.reset_index(drop=True, inplace=True)
        self.atom_type_set |= outer_sdat.atom_type_set
        for dim in range(3):
            self.cell[dim] = max(self.cell[dim], outer_sdat.cell[dim])

    def count_bonds(self, cut_off) -> dict:
        return bond_analysis.count_bonds(self, cut_off)

    def get_atom_idx_from_mol(self, cut_off, target_mol):
        return molecule_analysis.get_atom_idx_from_mol(self, cut_off, target_mol)

    def count_mols(self, cut_off, lower_mol_limit=1, upper_mol_limit=10) -> dict:
        return molecule_analysis.count_mols(self, cut_off, lower_mol_limit, upper_mol_limit)
