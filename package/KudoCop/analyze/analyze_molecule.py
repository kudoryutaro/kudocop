import sys
from collections import deque
from tqdm import tqdm, trange
import pandas as pd
import numpy as np


class AnalyzeMolecule():
    def __init__():
        pass

    def count_mols(self, cut_off=0.5, bond_type='dumpbond', lower_mol_limit=1, upper_mol_limit=10, condition=None) -> dict:
        """
        分子数をカウントする関数

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

        lower_mol_limit : int
            分子内の原子数がlower_mol_limit以上の分子のみをカウントする

        upper_mol_limit : int
            分子内の原子数がupper_mol_limit以下の分子のみをカウントする

        condition : function
            カウントしたい分子の条件を指定する関数

        Returns
        -------
        mol_counter : dict
            keyが分子、valueがその分子の個数となるdict
        """
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
        """
        分子を指定し、その分子のindexが入ったリストを返す関数
        例えば水分子の場合、
        [[315, 634, 123], [13, 75, 4124], [43, 1324, 63]]
        というリストが帰った場合、[315, 634, 123]が一つの水分子になっている。

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

        target_mol : tuple
            調べたい分子の中に何の原子が何個あるのかを指定するタプル。
            たとえば水分子を指定する時は
            {'C' : 1, 'H' : 2, 'O' : 3, 'N' : 4, 'Si' : 5}
            の場合は
            target_mol = (0, 0, 2, 1, 0, 0)
            となる。注意点はtarget_molのうちのtarget_mol[0]には0が入っていなければならない
            今後target_molは文字列で指定できるようにしたい。'H2O'みたいな

        condition : function
            カウントしたい分子の条件を指定する関数

        Returns
        -------
        atom_idx_from_mol : list
            それぞれにまとまった分子のindexが入ったリストのリスト
        """
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

    def get_connected_atoms_from_atom_idx(self, cut_off=0.5, bond_type='dumpbond', start_atom_idx=0, res_type='list') -> list:
        ''''
        start_atom_idxの原子につながっているすべての原子のインデックスが入ったリストを返す

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

        start_atom_idx : int
            つながっている原子すべてを数えるが、その始まりの原子のインデックス

        res_type : str
            res_type == 'list' の時は、つながっている原子すべてのインデックスが入ったリストを返す
            res_type == 'bool' の時は、つながっている原子はTrue,つながっていない原子は
            Falseとなっているnp_arrayを返す

        Returns
        -------
        connected_atoms_idxs : list or np.array
            res_type == 'list' の時は、つながっている原子すべてのインデックスが入ったリストを返す
            res_type == 'bool' の時は、つながっている原子はTrue,つながっていない原子は
            Falseとなっているnp_arrayを返す

        '''

        connect_list = self.get_connect_list(
            cut_off=cut_off, bond_type=bond_type)
        que = deque([start_atom_idx])

        if res_type == 'list':
            connected_atoms_idxs = set()
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
        elif res_type == 'bool':
            connected_atoms_bool = np.array([False] * self.get_total_atoms())
            while que:
                atom_idx = que.popleft()
                if connected_atoms_bool[atom_idx]:
                    continue
                connected_atoms_bool[atom_idx] = True
                for next_atom_idx in connect_list[atom_idx]:
                    if connected_atoms_bool[next_atom_idx]:
                        continue
                    que.append(next_atom_idx)
            return connected_atoms_bool

    def add_mol_q_to_atoms(self, cut_off=0.5, bond_type='dumpbond') -> None:
        """
        分子ごとの電荷を計算し、atomsにmol_qとして追加する

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
        """
        visited = [False] * self.get_total_atoms()
        sdat_q = self.atoms['q'].values
        sdat_mol_q = np.array([0.] * self.get_total_atoms())
        connect_list = self.get_connect_list(
            cut_off=cut_off, bond_type=bond_type)
        for atom_idx in range(self.get_total_atoms()):
            mol = []
            if visited[atom_idx]:
                continue
            que = deque()

            que.append(atom_idx)

            while que:
                current_atom_idx = que.popleft()
                if visited[current_atom_idx]:
                    continue
                visited[current_atom_idx] = True
                mol.append(current_atom_idx)
                for next_atom_idx in connect_list[current_atom_idx]:
                    if visited[next_atom_idx]:
                        continue
                    que.append(next_atom_idx)
            q = 0
            for atom_idx_in_mol in mol:
                q += sdat_q[atom_idx_in_mol]
            q = round(q, 6)
            for atom_idx_in_mol in mol:
                sdat_mol_q[atom_idx_in_mol] = q
        self.atoms['mol_q'] = sdat_mol_q


class AnalyzeMoleculeForSDats():
    def __init__():
        pass

    def count_mols(self, cut_off=0.5, bond_type='dumpbond', lower_mol_limit=1, upper_mol_limit=10, condition=None) -> pd.DataFrame:
        """
        分子数をカウントする関数

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

        lower_mol_limit : int
            分子内の原子数がlower_mol_limit以上の分子のみをカウントする

        upper_mol_limit : int
            分子内の原子数がupper_mol_limit以下の分子のみをカウントする

        condition : function
            カウントしたい分子の条件を指定する関数

        Returns
        -------
        df_mol_counter : DataFrame
            indexがstep_num、columnsが分子、dataがその分子の個数となるDataFrame
        """
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

    def add_mol_q_to_atoms(self, cut_off=0.5, bond_type='dumpbond') -> None:
        """
        分子ごとの電荷を計算し、atomsにmol_qとして追加する

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
        """
        connect_lists = self.get_connect_lists(
            cut_off=cut_off, bond_type=bond_type)
        for step_idx in range(len(self.step_nums)):

            visited = [False] * self.get_total_atoms()
            sdat_q = self.atoms[step_idx]['q'].values
            sdat_mol_q = np.array([0.] * self.get_total_atoms())
            for atom_idx in range(self.get_total_atoms()):
                mol = []
                if visited[atom_idx]:
                    continue
                que = deque()

                que.append(atom_idx)

                while que:
                    current_atom_idx = que.popleft()
                    if visited[current_atom_idx]:
                        continue
                    visited[current_atom_idx] = True
                    mol.append(current_atom_idx)
                    for next_atom_idx in connect_lists[step_idx][current_atom_idx]:
                        if visited[next_atom_idx]:
                            continue
                        que.append(next_atom_idx)
                q = 0
                for atom_idx_in_mol in mol:
                    q += sdat_q[atom_idx_in_mol]
                q = round(q, 6)
                for atom_idx_in_mol in mol:
                    sdat_mol_q[atom_idx_in_mol] = q
            self.atoms[step_idx]['mol_q'] = sdat_mol_q
