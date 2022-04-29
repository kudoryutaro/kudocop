from .atomic_info import atomic_weight
import pandas as pd
from .atomic_info import atomic_weight


class SimulationDat():
    def __init__(self):
        self.total_particle = 0
        self.particles = {}
        self.dcell = [0.0] * 3
        self.cell = [0.0] * 3
        self.newcell = [0.0] * 3
        self.elem_to_type = self._get_elem_to_type()
        self.type_to_elem = self._get_swap_dict(self.elem_to_type)
        self.type_set = set()
        self.type_to_mass = {_type:atomic_weight.get_weight(elem) for elem,_type in self.elem_to_type.items()}
        self.colum2suffix = {'pos':('x', 'y', 'z'), 'velo':('vx', 'vy', 'vz'), 'force':('fx', 'fy', 'fz')}
        self.connect_list = []

        #variables for config.rd
        self.mpigrid = [1] * 3
        self.ompgrid = [1] * 3
        self.cutoff = 2.5
        self.margin = 0.3
        self.total_step = 0
        self.file_step = 1000

        self.flagdecomp = False
        self.flagconnect = False

        #calclation infomation in input files
        self.fix_info = []
        self.move_info = []
        self.press_info = []
        self.sumforce_info = []
        self.thermofree_info= []
        self.wall_info = []

    @staticmethod
    def _get_elem_to_type():
        elem_to_type = {}
        for ind, key in enumerate(atomic_weight.atom_symbol_to_weight):
            elem_to_type[key] = ind + 1
        return elem_to_type

    @staticmethod
    def _get_swap_dict(d):
        return {v: k for k, v in d.items()}