import sys
from collections import deque
from tqdm import tqdm, trange
import pandas as pd
import numpy as np


class AnalyzeMolecule():
    def __init__():
        pass

    def count_mols(self, cut_off=0.5, bond_type='dumpbond', lower_mol_limit=1, upper_mol_limit=10, condition=None) -> dict:
        mol_counter = dict()
        connect_list = self.get_connect_list(
            cut_off=cut_off, bond_type=bond_type)

        visited = [0] * self.get_total_atoms()
        if condition is None:
            target_atoms = np.array([True] * self.get_total_atoms())
        else:
            target_atoms = condition(self)
        # number of types 種類数
        atom_type_set = self.get_atom_type_set()
        type_num_max = max(atom_type_set)
        sdat_atoms_type = self.atoms['type'].values

        for start_atom_idx in range(self.get_total_atoms()):
            if visited[start_atom_idx]:
                continue
            # 調べる範囲内に一つでも原子が入っている分子はすべてカウントする
            if not target_atoms[start_atom_idx]:
                continue
            current_mol_counter = [0] * (type_num_max + 1)

            # 幅優先探索
            que = deque([start_atom_idx])
            while que:
                current_atom_idx = que.popleft()
                if visited[current_atom_idx]:
                    continue
                visited[current_atom_idx] = 1
                current_mol_counter[sdat_atoms_type[current_atom_idx]] += 1
                for next_atom_idx in connect_list[current_atom_idx]:
                    if visited[next_atom_idx]:
                        continue
                    que.append(next_atom_idx)

            if lower_mol_limit <= sum(current_mol_counter) <= upper_mol_limit:
                current_mol_counter_tuple = tuple(current_mol_counter)
                if current_mol_counter_tuple not in mol_counter:
                    mol_counter[current_mol_counter_tuple] = 0
                mol_counter[current_mol_counter_tuple] += 1
        mol_counter_renamed = dict()
        for mol, count in mol_counter.items():
            mol_renamed = ''
            for atom_type, atom_count in enumerate(mol[1:], start=1):
                if atom_count != 0:
                    mol_renamed += f'{self.atom_type_to_symbol[atom_type]}{atom_count}'
            mol_counter_renamed[mol_renamed] = count

        return mol_counter_renamed

    def get_atom_idx_from_mol(self, cut_off=0.5, bond_type='dumpbond', target_mol=None, condition=None) -> list:
        if type(target_mol) != tuple:
            print('target_mol must be tuple')
            sys.exit(-1)
        if condition is None:
            target_atoms = np.array([True] * self.get_total_atoms())
        else:
            target_atoms = condition(self)
        atom_idx_from_mol = []
        visited = [0] * self.get_total_atoms()

        # number of types 種類数
        atom_type_set = self.get_atom_type_set()
        type_num_max = max(atom_type_set)
        if len(target_mol) != type_num_max + 1:
            print('target_mol\'s length does not match')
            sys.exit(-1)
        connect_list = self.get_connect_list(cut_off, bond_type)
        sdat_atoms_type = self.atoms['type'].values
        for start_atom_idx in range(self.get_total_atoms()):
            if visited[start_atom_idx]:
                continue
            # 調べる範囲内に一つでも原子が入っている分子はすべてカウントする
            if not target_atoms[start_atom_idx]:
                continue
            current_mol_counter = [0] * (type_num_max + 1)
            current_mol = []

            que = deque([start_atom_idx])

            while que:
                current_atom_idx = que.popleft()
                if visited[current_atom_idx]:
                    continue
                visited[current_atom_idx] = 1
                current_mol.append(current_atom_idx)
                current_mol_counter[sdat_atoms_type[current_atom_idx]] += 1

                for next_atom_idx in connect_list[current_atom_idx]:
                    if visited[next_atom_idx]:
                        continue
                    que.append(next_atom_idx)

            if tuple(current_mol_counter) == target_mol:
                atom_idx_from_mol.append(current_mol)

        return atom_idx_from_mol

    def get_connected_atoms_from_atom_idx(self, cut_off, start_atom_idx):
        ''''
        原子のインデックスがstart_atom_idxの原子につながっている原子のインデックスが入ったリストを返す
        '''

        connect_list = self.get_connect_list(cut_off)
        connected_atoms_idxs = set()
        que = deque([start_atom_idx])

        while que:
            atom_idx = que.popleft()
            if atom_idx in connected_atoms_idxs:
                continue
            connected_atoms_idxs.add(atom_idx)

            for next_atom_idx in connect_list[atom_idx]:
                if next_atom_idx in connected_atoms_idxs:
                    continue
                que.append(next_atom_idx)

        return list(connected_atoms_idxs)


class AnalyzeMoleculeForSDats():
    def __init__():
        pass

    def count_mols(self, cut_off=0.5, bond_type='dumpbond', lower_mol_limit=1, upper_mol_limit=10, condition=None) -> pd.DataFrame:
        dfs_count_mols = []
        atom_type_set = self.get_atom_type_set()
        type_num_max = max(atom_type_set)
        connect_lists = self.get_connect_lists(
            cut_off=cut_off, bond_type=bond_type)
        for step_idx, step_num in enumerate(tqdm(self.step_nums, desc='[counting mols]')):
            mol_counter = dict()
            visited = [0] * self.get_total_atoms()
            if condition is None:
                target_atoms = np.array([True] * self.get_total_atoms())
            else:
                target_atoms = condition(self, step_idx)
            # number of types 種類数

            sdat_atoms_type = self.atoms[step_idx]['type'].values

            for start_atom_idx in range(self.get_total_atoms()):
                if visited[start_atom_idx]:
                    continue
                if not target_atoms[start_atom_idx]:
                    continue
                current_mol_counter = [0] * (type_num_max + 1)

                # 幅優先探索
                que = deque([start_atom_idx])
                while que:
                    current_atom_idx = que.popleft()
                    if visited[current_atom_idx]:
                        continue
                    visited[current_atom_idx] = 1
                    current_mol_counter[sdat_atoms_type[current_atom_idx]] += 1
                    for next_atom_idx in connect_lists[step_idx][current_atom_idx]:
                        if visited[next_atom_idx]:
                            continue
                        que.append(next_atom_idx)

                if lower_mol_limit <= sum(current_mol_counter) <= upper_mol_limit:
                    current_mol_counter_tuple = tuple(current_mol_counter)
                    if current_mol_counter_tuple not in mol_counter:
                        mol_counter[current_mol_counter_tuple] = 0
                    mol_counter[current_mol_counter_tuple] += 1

            dfs_count_mols.append(pd.DataFrame(mol_counter, [step_num]))
        df_count_mols = pd.concat(dfs_count_mols).fillna(0)
        df_count_mols.columns = list(df_count_mols.columns)

        def change_column_name(column):
            column_res = []
            for atom_type, atom_num in enumerate(column[1:], start=1):
                if atom_num == 0:
                    continue
                if atom_type not in self.atom_type_to_symbol:
                    column_res.append(f'X{atom_num}')
                else:
                    column_res.append(
                        f'{self.atom_type_to_symbol[atom_type]}{atom_num}')
            return ''.join(column_res)
        df_count_mols = df_count_mols.rename(columns=change_column_name)

        return df_count_mols
