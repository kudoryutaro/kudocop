import pandas as pd
from .io import read_para
from .io import read_input
from .io import read_config
from .io import read_dumppos
from .io import read_dumpbond
from .io import read_xyz


class SimulationDat():
    def __init__(self):
        self.atoms = None
        self.cell = [None] * 3
        self.atom_symbol_to_type = None
        self.atom_type_to_symbol = None
        self.atom_type_to_mass = None
        self.atom_type_set = set()
        self.bondorder_list = None
        self.connect_list = None

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

    # METHODS

    def get_total_atoms(self):
        return len(self.atoms)

    def set_decomp(self):
        clx = float(1.0 / (self.ompgrid[0] * self.mpigrid[0]))
        cly = float(1.0 / (self.ompgrid[1] * self.mpigrid[1]))
        clz = float(1.0 / (self.ompgrid[2] * self.mpigrid[2]))
        self.decompx = [clx] * self.ompgrid[0] * self.mpigrid[0]
        self.decompy = [cly] * self.ompgrid[1] * self.mpigrid[1]
        self.decompz = [clz] * self.ompgrid[2] * self.mpigrid[2]
