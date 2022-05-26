import sys
import pandas as pd
import numpy as np
from .io.export_file import ExportFile
from .io.import_file import ImportFile
from .analyze.analyze import Analyze
from .analyze import neighbor


class SimulationDat(
    ImportFile,
    ExportFile,
    Analyze
):
    def __init__(self):
        self.atoms = None
        self.cell = [None] * 3
        self.bondorder_list = None
        self.bondorder_connect_list = None
        self.connect_list = None
        self.connect_list_cut_off = None
        self.connect_list_from_dumppos = None
        self.connect_list_cut_off_from_dumppos = None

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

    # METHODS\
    def __getitem__(self, key):
        return self.atoms[key]

    def __setitem__(self, key, val):
        self.atoms[key] = val

    def __len__(self):
        return self.get_total_atoms()

    def get_total_atoms(self):
        if self.atoms is not None:
            return len(self.atoms)
        elif self.bondorder_list is not None:
            return len(self.bondorder_list)
        else:
            print('Import file first')
            sys.exit(-1)

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

    def get_connect_list_from_dumppos(self, cut_off):
        if self.connect_list_cut_off_from_dumppos != cut_off or self.connect_list_from_dumppos is None:
            self.__create_connect_list_from_dumppos(cut_off)
        self.connect_list_cut_off_from_dumppos = cut_off
        return self.connect_list_from_dumppos

    def __create_connect_list_from_dumppos(self, cut_off):
        self.particles = dict()
        self.particles['pos'] = self.atoms[['x', 'y', 'z']].values
        self.total_particle = self.get_total_atoms()
        self.newcell = self.cell
        self.connect_list_from_dumppos = neighbor.make_neighbor(self, cut_off)
        return self.connect_list_from_dumppos

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
