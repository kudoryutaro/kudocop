import sys
from collections import deque
from tqdm import tqdm, trange
import pandas as pd


class AnalyzeMolecule():
    def __init__():
        pass

    def count_mols(self, cut_off, lower_mol_limit=1, upper_mol_limit=10) -> dict:
        mol_counter = dict()
        self.get_connect_list(cut_off)

        visited = [0] * self.get_total_atoms()

        # number of types 種類数
        type_num_max = max(self.atom_type_set)

        self_atoms_type = self.atoms['type'].values

        for start_atom_idx in range(self.get_total_atoms()):
            if visited[start_atom_idx]:
                continue

            current_mol_counter = [0] * (type_num_max + 1)

            # 幅優先探索
            que = deque([start_atom_idx])
            while que:
                current_atom_idx = que.popleft()
                if visited[current_atom_idx]:
                    continue
                visited[current_atom_idx] = 1
                current_mol_counter[self_atoms_type[current_atom_idx]] += 1
                for next_atom_idx in self.connect_list[current_atom_idx]:
                    if visited[next_atom_idx]:
                        continue
                    que.append(next_atom_idx)

            if lower_mol_limit <= sum(current_mol_counter) <= upper_mol_limit:
                current_mol_counter_tuple = tuple(current_mol_counter)
                if current_mol_counter_tuple not in mol_counter:
                    mol_counter[current_mol_counter_tuple] = 0
                mol_counter[current_mol_counter_tuple] += 1

        return mol_counter

    # def count_mols(selfs, cut_off, lower_mol_limit, upper_mol_limit, rename_columns=True) -> pd.DataFrame:
    #     dfs_count_mols = []
    #     for step_idx, step_num in enumerate(tqdm(selfs.step_nums)):
    #         counter = selfs.data[step_idx].count_mols(
    #             cut_off, lower_mol_limit, upper_mol_limit)
    #         dfs_count_mols.append(pd.DataFrame(counter, [step_num]))
    #     df_count_mols = pd.concat(dfs_count_mols).fillna(0)
    #     df_count_mols.columns = list(df_count_mols.columns)

    #     def change_column_name(column):
    #         column_res = []
    #         for atom_type, atom_num in enumerate(column[1:], start=1):
    #             if atom_num == 0:
    #                 continue
    #             if atom_type not in selfs.data[0].atom_type_to_symbol:
    #                 column_res.append(f'X{atom_num}')
    #             else:
    #                 column_res.append(
    #                     f'{selfs.data[0].atom_type_to_symbol[atom_type]}{atom_num}')
    #         return ''.join(column_res)
    #     if rename_columns:
    #         df_count_mols = df_count_mols.rename(columns=change_column_name)

    #     return df_count_mols

    def get_atom_idx_from_mol(self, cut_off, target_mol):
        if type(target_mol) != tuple:
            print('target_mol must be tuple')
            sys.exit(-1)

        atom_idx_from_mol = []
        self.get_connect_list(cut_off)
        visited = [0] * self.get_total_atoms()

        # number of types 種類数
        type_num_max = max(self.atom_type_set)
        if len(target_mol) != type_num_max + 1:
            print('target_mol\'s length does not match')
            sys.exit(-1)

        self_atoms_type = self.atoms['type'].values
        for start_atom_idx in range(self.get_total_atoms()):
            if visited[start_atom_idx]:
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
                current_mol_counter[self_atoms_type[current_atom_idx]] += 1

                for next_atom_idx in self.connect_list[current_atom_idx]:
                    if visited[next_atom_idx]:
                        continue
                    que.append(next_atom_idx)

            if tuple(current_mol_counter) == target_mol:
                atom_idx_from_mol.append(current_mol)

        return atom_idx_from_mol
