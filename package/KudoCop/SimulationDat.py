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

        self.connect_list_cut_off_from_dumpbond = None
        self.connect_list_from_dumpbond = None

        self.connect_list_from_dumppos = None
        self.connect_list_cut_off_from_dumppos = None

        self.connect_list_from_dumpbond_cg = None

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

        self.sumforce = None

    # METHODS\
    def __getitem__(self, key) -> pd.DataFrame:
        """
        sdat.atoms[column]をsdat[column]と省略して書くことが出来る。
        """
        return self.atoms[key]

    def __setitem__(self, key, val) -> None:
        self.atoms[key] = val

    def __len__(self) -> int:
        """
        sdat.get_total_atoms()をlen(sdat)と省略して書くことが出来る。
        """
        return self.get_total_atoms()

    def get_total_atoms(self) -> int:
        """
        全原子数を返す関数。
        """
        if self.atoms is not None:
            return len(self.atoms)
        elif self.bondorder_list is not None:
            return len(self.bondorder_list)
        elif self.connect_list_from_dumpbond_cg is not None:
            return len(self.connect_list_from_dumpbond_cg)
        else:
            print('Import file first')
            sys.exit(-1)

    def get_atom_type_set(self) -> set:
        """
        系内の原子のtypeのsetを返す関数
        """
        return set(self.atoms['type'])

    def wrap_particles(self) -> None:
        """
        セルの外にはみ出している原子をセルの中に入れる。
        """

        if self.cell is None:
            print('set sdat.cell')
            sys.exit(-1)
        if 0 in self.cell:
            print('cell size must not be 0')
            sys.exit(-1)
        self.atoms[['x', 'y', 'z']] %= self.cell

    def get_connect_list(self, cut_off=0.5, bond_type='dumpbond') -> list:
        """
        connect_listを返す関数。

        Parameters
        ----------
        cut_off : float
            bond_type == 'dumppos' の時はcut_offの単位はÅ
            ある原子からcut_off以下の距離にある原子は結合しているとみなす

            bond_type == 'dumpbond' の時はcut_offの単位はbond order
            ある原子とある原子のbond orderの和がcut_off以上のときに結合しているとみなす

            bond_type == 'dumpbond_cg' の時はcut_offは不必要

        bond_type : str
            bond_type == 'dumppos' の時はconnect_listはdumpposから生成される
            bond_type == 'dumpbond' の時connect_listはdumpbondから生成される
            bond_type == 'dumpbond_cg' の時connect_listはdumpbond_cgから生成される

        Returns
        -------
        connect_list : list
            原子の隣接リスト
            connect_list[atom_idx] -> [neibour_atom_idx1, neibour_atom_idx2, ...]

        """
        if bond_type == 'dumppos':
            return self.get_connect_list_from_dumppos(cut_off)
        elif bond_type == 'dumpbond':
            return self.get_connect_list_from_dumpbond(cut_off)
        elif bond_type == 'dumpbond_cg':
            return self.connect_list_from_dumpbond_cg
        else:
            print('unsupported bond_type')
            sys.exit(-1)

    def get_connect_list_from_dumpbond(self, cut_off) -> list:
        """
        dumpbondからconnect_listを生成する関数。
        以前にcut_offが同じconnect_listが生成されていたら新しく生成せずに前に生成したconnect_listを返す

        Parameters
        ----------
        cut_off : float
            cut_offの単位はbond order
            ある原子とある原子のbond orderの和がcut_off以上のときに結合しているとみなす

        Returns
        -------
        connect_list : list
            原子の隣接リスト
            connect_list[atom_idx] -> [neibour_atom_idx1, neibour_atom_idx2, ...]

        """
        if self.connect_list_cut_off_from_dumpbond != cut_off or self.connect_list is None:
            self.__create_connect_list_from_dumpbond(cut_off)
        self.connect_list_cut_off_from_dumpbond = cut_off
        return self.connect_list

    def __create_connect_list_from_dumpbond(self, cut_off) -> None:
        """
        dumpbondからconnect_listを生成する関数。

        Parameters
        ----------
        cut_off : float
            cut_offの単位はbond order
            ある原子とある原子のbond orderの和がcut_off以上のときに結合しているとみなす
        """
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

    def get_connect_list_from_dumppos(self, cut_off) -> list:
        """
        dumpposからconnect_listを生成する関数。
        以前にcut_offが同じconnect_listが生成されていたら新しく生成せずに前に生成したconnect_listを返す

        Parameters
        ----------
        cut_off : float
            cut_offの単位はÅ
            ある原子からcut_off以下の距離にある原子は結合しているとみなす

        Returns
        -------
        connect_list : list
            原子の隣接リスト
            connect_list[atom_idx] -> [neibour_atom_idx1, neibour_atom_idx2, ...]

        """
        if self.connect_list_cut_off_from_dumppos != cut_off or self.connect_list_from_dumppos is None:
            self.__create_connect_list_from_dumppos(cut_off)
        self.connect_list_cut_off_from_dumppos = cut_off
        return self.connect_list_from_dumppos

    def __create_connect_list_from_dumppos(self, cut_off) -> list:
        """
        dumpposからconnect_listを生成する関数。

        Parameters
        ----------
        cut_off : float
            cut_offの単位はÅ
            ある原子からcut_off以下の距離にある原子は結合しているとみなす
        """
        self.particles = dict()
        self.particles['pos'] = self.atoms[['x', 'y', 'z']].values
        self.total_particle = self.get_total_atoms()
        self.newcell = self.cell
        self.connect_list_from_dumppos = neighbor.make_neighbor(self, cut_off)
        return self.connect_list_from_dumppos

    def replicate_atoms(self, replicate_directions=[1, 1, 1]) -> None:
        """
        x, y, z 方向にセルを複製する関数

        Parameters
        ----------
        replicate_directions : list
            x, y, z方向に何倍するかを指定する。
            例えばx方向に2倍,y方向に3倍,z方向に4倍したい時は
            replicate_directions = [2, 3, 4] とする

        """
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

    def concat_atoms(self, outer_sdat) -> None:
        """
        sdatとouter_sdatを結合する関数

        Parameters
        ----------
        outer_sdat : SimulationDat
            取り入れたいsdatを指定する。

        """
        self.atoms = pd.concat([self.atoms, outer_sdat.atoms])
        self.atoms.reset_index(drop=True, inplace=True)
        for dim in range(3):
            self.cell[dim] = max(self.cell[dim], outer_sdat.cell[dim])

    def delete_atoms(self, condition, reindex):
        """
        条件に当てはまる原子を削除する

        Parameters
        ----------
        condition : function
            削除したい原子の条件を指定する関数
        reindex : bool
            reindex == Trueの時は原子のid(idx)は新しく割り振られる
            reindex == Falseの時は原子のid(idx)は新しく割り振られず、削除前のものと変わらない
        """
        if callable(condition):
            target_atoms = condition(self)
            self.atoms = self.atoms[~target_atoms]
        else:
            self.atoms = self.atoms[~condition]
        if reindex:
            self.atoms.reset_index(drop=True, inplace=True)
