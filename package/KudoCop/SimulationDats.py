import sys
import pandas as pd
import numpy as np
import os
from tqdm import tqdm, trange
from .SimulationDat import SimulationDat
from .io.import_dumppos import ImportDumppos
from .io.import_dumpbond import ImportDumpbond
from .io.import_para import ImportPara
from .analysis import bond_analysis
from .analysis import atom_analysis
from .analysis import molecule_analysis


class SimulationDats():
    def __init__(self, para_file_name: str, dir_name=None, step_nums=None):
        if dir_name is None:
            self.dir_name = os.getcwd()
        else:
            self.dir_name = dir_name
        self.para_file_name = para_file_name
        file_names_in_current_dir = os.listdir(dir_name)
        if step_nums is None:
            self.step_nums = []
            for file_name in file_names_in_current_dir:
                if len(file_name) >= 9 and file_name[:9] == 'dump.pos.':
                    self.step_nums.append(int(file_name[9:]))
        else:
            self.step_nums = step_nums
        self.step_nums.sort()

        self.step_num_to_step_idx = {
            step_num: step_idx for step_idx, step_num in enumerate(self.step_nums)
        }

        self.data = [SimulationDat() for _ in range(len(self.step_nums))]

        for step_idx, step_num in enumerate(tqdm(self.step_nums)):
            dumppos_importer = ImportDumppos()
            dumppos_importer.import_file(
                self.data[step_idx], f'{self.dir_name}/dump.pos.{step_num}')
            dumpbond_importer = ImportDumpbond()
            dumpbond_importer.import_file(
                self.data[step_idx], f'{self.dir_name}/dump.bond.{step_num}')
            para_importer = ImportPara()
            para_importer.import_file(
                self.data[step_idx], f'{self.dir_name}/{self.para_file_name}')

    def count_mols(self, cut_off, lower_mol_limit=1, upper_mol_limit=10, rename_columns=True) -> pd.DataFrame:
        return molecule_analysis.count_mols_sdats(self, cut_off, lower_mol_limit, upper_mol_limit, rename_columns)

    def count_bonds(self, cut_off):
        return bond_analysis.count_bonds_sdats(self, cut_off)

    def to_dumppos(self, output_folder, out_columns=None):
        try:
            os.makedirs(output_folder)
        except FileExistsError:
            pass

        for step_idx, step_num in tqdm(enumerate(self.step_nums)):
            self.data[step_idx].to_dumppos(
                f'{output_folder}/dump.pos.{step_num}', time_step=step_num, out_columns=out_columns)
