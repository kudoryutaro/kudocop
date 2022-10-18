import sys
import pandas as pd
import numpy as np
import os
from tqdm import tqdm, trange
from .io.import_para import ImportPara
from .io.import_dumpposes import ImportDumpposes
from .io.import_dumpbonds_cg import ImportDumpbondsCG
from .io.import_dumpbonds import ImportDumpbonds
from .io.export_dumpposes import ExportDumpposes
from .analyze.analyze import AnalyzeForSDats
from .SimulationDat import SimulationDat
from .analyze import neighbor
from .io.export_dp_system_multi_frames import ExportDPSystemMultiFrames
from .analyze.dmol3 import DMol3KudoCopForSDats
from .io.import_outmol import ImportOutmolForSDats
from .io.import_lammps_dumppos import ImportLammpsDumppos

class SimulationDats(
    ImportPara,
    ImportDumpposes,
    ImportDumpbonds,
    ImportDumpbondsCG,
    ImportOutmolForSDats,
    ExportDumpposes,
    AnalyzeForSDats,
    DMol3KudoCopForSDats,
    ExportDPSystemMultiFrames,
    ImportLammpsDumppos,
):
    """シミュレーションしたデータを読み込み、書き込み、分析するためのクラス
    一度に複数のdump.posとdump.bondを読み込む

    Attributes
    ----------
    atoms : list
        atoms[step_idx] : pd.DataFrame
        原子の座標, type, 電荷などを含むpandasのDataFrame
    cell : list
        cellの大きさが入ったlist
        cell[0]:x方向, cell[1]:y方向, cell[2]:z方向
    bondorder_lists : list
        dump.bondを読み込んだbond orderのリスト
    bondorder_connect_lists : list
        dump.bondを読み込んだconnect_listのリスト
    connect_list_cut_off_from_dumpbonds : float
        一度connect_listを作成した時のcut_off
        同じcut_offで何度も同じconnect_listを作らないように保持する
    connect_list_cut_off_from_dumpposes : float
        一度connect_listを作成した時のcut_off
        同じcut_offで何度も同じconnect_listを作らないように保持する
    atom_symbol_to_type : dict
        原子のシンボルをkey, 原子のtypeをvalueとするdict
    atom_type_to_symbol : dict
        原子のtypeをkey, 原子のシンボルをvalueとするdict
    atom_type_to_mass : dict
        原子のtypeをkey, 原子の質量(g/mol)をvalueとするdict
    
    """
    def __init__(self, para_file_name='para.rd', dir_name=None, import_dumpposes_flag=True, import_dumpbonds_flag=True, step_nums=None, skip_num=None,bond_type='dumpbond',max_coordination_num=6):

        self.cell = [None] * 3
        self.bondorder_lists = None
        self.bondorder_connect_lists = None

        self.connect_lists_from_dumpbonds = None
        self.connect_list_cut_off_from_dumpbonds = None

        self.connect_lists_from_dumpposes = None
        self.connect_list_cut_off_from_dumpposes = None

        self.potential_energy = []

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
            if bond_type == 'dumpbond':
                self.import_dumpbonds()
            elif bond_type == 'dumpbond_cg':
                self.import_dumpbonds_cg(max_coordination_num)

    def __getitem__(self, key):
        return self.atoms[key]

    def __len__(self):
        return len(self.step_nums)

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

    def get_connect_lists(self,cut_off=0.5, bond_type='dumpbond'):
        if bond_type == 'dumppos':
            return self.get_connect_lists_from_dumpposes(cut_off)
        elif bond_type == 'dumpbond':
            return self.get_connect_lists_from_dumpbonds(cut_off)
        elif bond_type == 'dumpbond_cg':
            return self.get_connect_lists_from_dumpbonds()
        else:
            print('unsupported bond_type')
            sys.exit(-1)
        
    def get_connect_lists_from_dumpposes(self,cut_off):
        if self.connect_list_cut_off_from_dumpposes != cut_off or self.connect_lists_from_dumpposes is None:
            self.__create_connect_lists_from_dumpposes(cut_off)
        self.connect_list_cut_off_from_dumpposes = cut_off
        return self.connect_lists_from_dumpposes
    
    def __create_connect_lists_from_dumpposes(self,cut_off):
        if cut_off is None:
            print('cut_off is not defined')
            sys.exit(-1)

        self.connect_lists_from_dumpposes = [None for _ in range(len(self.step_nums))]

        for step_idx in trange(len(self.step_nums), desc='[creating connect_lists]'):
            sdat = SimulationDat()
            sdat.particles = dict()
            sdat.particles['pos'] = self.atoms[step_idx][['x', 'y', 'z']].values
            sdat.total_particle = self.get_total_atoms()
            sdat.newcell = self.cell
            self.connect_lists_from_dumpposes[step_idx] = neighbor.make_neighbor(sdat, cut_off)
        return self.connect_lists_from_dumpposes


    def get_connect_lists_from_dumpbonds(self, cut_off):
        if self.connect_list_cut_off_from_dumpbonds != cut_off or self.connect_lists_from_dumpbonds is None:
            self.__create_connect_lists_from_dumpbonds(cut_off)
        self.connect_list_cut_off_from_dumpbonds = cut_off
        return self.connect_lists_from_dumpbonds

    def __create_connect_lists_from_dumpbonds(self, cut_off):
        if cut_off is None:
            print('cut_off is not defined')
            sys.exit(-1)
        if self.bondorder_lists is None:
            print('bondorder_list is not defined')
            print('Import dumpbond first')
            sys.exit(-1)
        self.connect_lists_from_dumpbonds = [[[] for _ in range(
            self.get_total_atoms())] for _ in range(len(self.step_nums))]

        for step_idx in trange(len(self.step_nums), desc='[creating connect_lists]'):
            for atom_idx, (neibour_idxs, bondorder_list) in enumerate(zip(self.bondorder_connect_lists[step_idx], self.bondorder_lists[step_idx])):
                for neibour_idx, bond_l in zip(neibour_idxs, bondorder_list):
                    if bond_l[-1] >= cut_off:
                        self.connect_lists_from_dumpbonds[step_idx][atom_idx].append(
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


    def reshape_bondorder_lists_cutoff(self, cut_off):
        """
        dumpbondsから作成したbondorder_lists内の、bond orderの和がcut_off未満の行を削除する関数
        
        Parameters
        ----------
        cut_off : float
            cut_offの単位はbond order
            ある原子とある原子のbond orderの和がcut_off未満のときの結合を除外する
        judge_cutoff : list
            ある原子([atom_idx])に着目した場合、bondorder_lists[atom_idx]の中で
            bond orderの和がcut_off未満である結合の、bondorder_lists[atom_idx]内での順番を保存する配列。

        """
        if cut_off is None:
            print('cut_off is not defined')
            sys.exit(-1)
        if self.bondorder_lists is None:
            print('bondorder_list is not defined')
            print('Import dumpbond first')
            sys.exit(-1)
        for step_idx in range(len(self.step_nums)):
            for atom_idx in range(len(self.atoms[0])):
                judge_cutoff = []
                for connected_atom_row in range(len(self.bondorder_lists[step_idx][atom_idx])):
                    if self.bondorder_lists[step_idx][atom_idx][connected_atom_row][-1] < cut_off:
                        judge_cutoff.append(connected_atom_row)
                self.bondorder_lists[step_idx][atom_idx] = np.delete(self.bondorder_lists[step_idx][atom_idx], judge_cutoff, 0)
        return self.bondorder_lists

    def concat_sdats(self, sdats_list:list):
        """sdatsを結合する
        Parameters
        ----------
            sdats_list : list of sdats
                結合するsdatsのリスト, 
        Note
        ----
            concat_sdatsメソッドを使用するsdatsは空のを使う
        """

        for outer_sdats in sdats_list:
            if self.cell[0] is None and self.cell[1] is None and self.cell[2] is None:
                if outer_sdats.cell[0] is not None and outer_sdats.cell[1] is not None and outer_sdats.cell[2] is not None:
                    self.cell[0] = outer_sdats.cell[0]
                    self.cell[1] = outer_sdats.cell[1]
                    self.cell[2] = outer_sdats.cell[2]

            for step_idx in range(len(outer_sdats.step_nums)):
                self.atoms.append(outer_sdats.atoms[step_idx])
                self.potential_energy.append(outer_sdats.potential_energy[step_idx])
        
        self.step_nums = list(range(len(self.atoms)))
        self.step_num_to_step_idx = {
           step_num:step_idx for step_idx, step_num in enumerate(range(len(self.step_nums)))
        }

    def add_mass_to_atoms(self):
        """原子の質量の列をatomsに追加する
        単位はg/mol

        Parameters
        ----------
            None
        """
        for step_idx in range(len(self.step_nums)):
            self.atoms[step_idx]['mass'] = self.atoms[step_idx]['type'].map(self.atom_type_to_mass)

    def add_force_to_atoms(self):
        """それぞれの原子の加速度からそれぞれの原子に働く力を求めてatomsに追加する
        単位はeV/Å
        Parameters
        ----------
            None
        """
        self.add_mass_to_atoms()
        hartree = 27.211386024367243
        bohr = 0.5291772105638411
        for step_idx in range(len(self.step_nums)):
            assert 'ax' in self.atoms[step_idx] and 'ay' in self.atoms[step_idx] and 'az' in self.atoms[step_idx], \
                'atomsに加速度がありません'
            
            # 単位換算 要検討
            for direction in ['x', 'y', 'z']:
                self.atoms[step_idx]['f' + direction] = \
                    self.atoms[step_idx]['a' + direction] * self.atoms[step_idx]['mass'] * hartree / bohr
