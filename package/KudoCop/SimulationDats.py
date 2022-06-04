import sys
import pandas as pd
import numpy as np
import os
from tqdm import tqdm, trange
from .io.import_para import ImportPara
from .io.import_dumpposes import ImportDumpposes
from .io.import_dumpbonds import ImportDumpbonds
from .io.export_dumpposes import ExportDumpposes
from .analyze.analyze import AnalyzeForSDats


class SimulationDats(
    ImportPara,
    ImportDumpposes,
    ImportDumpbonds,
    ExportDumpposes,
    AnalyzeForSDats
):
    def __init__(self, para_file_name='para.rd', dir_name=None, import_dumpposes_flag=True, import_dumpbonds_flag=True, step_nums=None, skip_num=None):
        self.cell = [None] * 3
        self.bondorder_lists = None
        self.bondorder_connect_lists = None
        self.connect_lists = None
        self.connect_list_cut_off = None



        # variables for para.rd
        self.atom_symbol_to_type = None
        self.atom_type_to_symbol = None
        self.atom_type_to_mass = None

        if dir_name is None:
            self.dir_name = os.getcwd()
        else:
            self.dir_name = dir_name
        file_names_in_current_dir = os.listdir(dir_name)
        if step_nums is None:
            self.step_nums = []
            for file_name in file_names_in_current_dir:
                if len(file_name) >= 9 and file_name[:9] == 'dump.pos.':
                    self.step_nums.append(int(file_name[9:]))

        else:
            self.step_nums = []
            for step_num in step_nums:
                self.step_nums.append(step_num)
        self.step_nums.sort()
        if skip_num is not None:
            self.step_nums = self.step_nums[::skip_num]

        self.step_num_to_step_idx = {
            step_num: step_idx for step_idx, step_num in enumerate(self.step_nums)
        }

        self.atoms = [None for _ in range(len(self.step_nums))]

        self.import_para(para_file_name)

        if import_dumpposes_flag:
            self.import_dumpposes()
        if import_dumpbonds_flag:
            self.import_dumpbonds()

    def get_total_atoms(self):
        if self.atoms[0] is not None:
            return len(self.atoms[0])
        elif self.bondorder_lists[0] is not None:
            return len(self.bondorder_lists[0])
        else:
            print('Import file first')
            sys.exit(-1)

    def get_atom_type_set(self):
        return set(self.atoms[0]['type'])

    def get_connect_lists(self, cut_off):
        if self.connect_list_cut_off != cut_off or self.connect_lists is None:
            self.__create_connect_lists(cut_off)
        self.connect_list_cut_off = cut_off
        return self.connect_lists

    def __create_connect_lists(self, cut_off):
        if cut_off is None:
            print('cut_off is not defined')
            sys.exit(-1)
        if self.bondorder_lists is None:
            print('bondorder_list is not defined')
            print('Import dumppos first')
            sys.exit(-1)
        self.connect_lists = [[[] for _ in range(
            self.get_total_atoms())] for _ in range(len(self.step_nums))]

        for step_idx in trange(len(self.step_nums), desc='[creating connect_lists]'):
            for atom_idx, (neibour_idxs, bondorder_list) in enumerate(zip(self.bondorder_connect_lists[step_idx], self.bondorder_lists[step_idx])):
                for neibour_idx, bond_l in zip(neibour_idxs, bondorder_list):
                    if bond_l[-1] >= cut_off:
                        self.connect_lists[step_idx][atom_idx].append(
                            neibour_idx)

    def delete_atoms(self, condition, reindex):
        # Trueの原子を削除する
        for step_idx in range(len(self.step_nums)):
            if callable(condition):
                target_atoms = condition(self, step_idx)
                self.atoms[step_idx] = self.atoms[step_idx][~target_atoms]
            else:
                self.atoms[step_idx] = self.atoms[step_idx][~condition]
            if reindex:
                self.atoms[step_idx].reset_index(drop=True, inplace=True)

    