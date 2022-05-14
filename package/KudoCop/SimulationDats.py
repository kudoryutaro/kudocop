import sys
import pandas as pd
import numpy as np
import os
from tqdm import tqdm, trange
from .io.import_para import ImportPara
from .io.import_dumpposes import ImportDumpposes
from .io.import_dumpbonds import ImportDumpbonds
from .analyze.analyze import AnalyzeForSDats


class SimulationDats(
    ImportPara,
    ImportDumpposes,
    ImportDumpbonds,
    AnalyzeForSDats
):
    def __init__(self, para_file_name='para.rd', dir_name=None, import_dumpposes_flag=True, import_dumpbonds_flag=True, step_nums=None):
        self.cell = [None] * 3
        self.bondorder_lists = None
        self.bondorder_connect_lists = None
        self.connect_lists = None
        self.connect_list_cut_off = None

        # atom set which exist in the system
        self.atom_type_set = set()

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

        for step_idx in trange(len(self.step_nums)):
            for atom_idx, (neibour_idxs, bondorder_list) in enumerate(zip(self.bondorder_connect_lists[step_idx], self.bondorder_lists[step_idx])):
                for neibour_idx, bond_l in zip(neibour_idxs, bondorder_list):
                    if bond_l[-1] >= cut_off:
                        self.connect_lists[step_idx][atom_idx].append(
                            neibour_idx)

    # def count_mols(self, cut_off, lower_mol_limit=1, upper_mol_limit=10, rename_columns=True) -> pd.DataFrame:
    #     return molecule_analysis.count_mols_sdats(self, cut_off, lower_mol_limit, upper_mol_limit, rename_columns)

    # def count_bonds(self, cut_off):
    #     return bond_analysis.count_bonds_sdats(self, cut_off)

    # def export_dumppos(self, output_folder, out_columns=None):
    #     try:
    #         os.makedirs(output_folder)
    #     except FileExistsError:
    #         pass

    #     for step_idx, step_num in tqdm(enumerate(self.step_nums)):
    #         self.data[step_idx].export_dumppos(
    #             f'{output_folder}/dump.pos.{step_num}', time_step=step_num, out_columns=out_columns)
