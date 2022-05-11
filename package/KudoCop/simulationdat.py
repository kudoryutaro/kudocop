import sys
import pandas as pd
import numpy as np
from .io import read_para
from .io import read_input
from .io import read_config
from .io import read_dumppos
from .io import read_dumpbond
from .io import read_xyz
from .io import to_input
from .io import to_dumppos
from .io import to_xyz


class SimulationDat():
    def __init__(self):
        self.atoms = None
        self.cell = [None] * 3
        self.bondorder_list = None
        self.connect_list = None

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

    # READ
    def read_input(self, ifn: str) -> None:
        reader = read_input.ReadInput()
        reader.read_file(self, ifn)

    def read_config(self, ifn: str) -> None:
        reader = read_config.ReadConfig()
        reader.read_file(self, ifn)

    def read_dumppos(self, ifn: str) -> None:
        reader = read_dumppos.ReadDumppos()
        reader.read_file(self, ifn)

    def read_dumpbond(self, ifn: str) -> None:
        reader = read_dumpbond.ReadDumpbond()
        reader.read_file(self, ifn)

    def read_xyz(self, ifn: str) -> None:
        reader = read_xyz.ReadXyz()
        reader.read_file(self, ifn)

    def read_para(self, ifn: str) -> None:
        reader = read_para.ReadPara()
        reader.read_file(self, ifn)

    # WEITE
    def to_input(self, ofn: str) -> None:
        writer = to_input.ToInput()
        writer.to_file(self, ofn)

    def to_dumppos(self, ofn: str, time_step=None, out_columns=None) -> None:
        writer = to_dumppos.ToDumppos()
        writer.to_file(self, ofn, time_step, out_columns)

    def to_xyz(self, ofn: str, out_columns=None, structure_name=None) -> None:
        writer = to_xyz.ToXyz()
        writer.to_file(self, ofn, out_columns, structure_name)

    # METHODS\
    def get_total_atoms(self):
        if self.atoms is not None:
            return len(self.atoms)
        elif self.bondorder_list is not None:
            return len(self.bondorder_list)
        else:
            print('Read file first')
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

    def __create_connect_list(self, cut_off):
        if cut_off is None:
            print('cut_off is not defined')
            sys.exit(-1)
        if self.bondorder_list is None:
            print('bondorder_list is not defined')
            print('Read dumppos first')
            sys.exit(-1)
        self.connect_list = [[] for _ in range(self.get_total_atoms())]
        for atom_idx, bondorder_list in enumerate(self.bondorder_list):
            for bond_l in bondorder_list:
                if bond_l[-1] >= cut_off:
                    self.connect_list[atom_idx].append(bond_l[0])

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
