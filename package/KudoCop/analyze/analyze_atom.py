from tqdm import tqdm, trange
import itertools
import pandas as pd
import numpy as np


class AnalyzeAtom():
    def __init__():
        pass

    def count_triplets(self, cut_off=0.5, bond_type='dumpbond', condition=None) -> dict:
        """
        ３体間の原子をカウントする関数
        例えば、水分子一つにH-O-Hは1個含まれている。
        メタン分子一つにはcombination(4,2) = 6より、H-C-Hは6個含まれている。

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

        condition : function
            3体間のうち、カウントしたい部分を指定する関数

        Returns
        -------
        count_triplets : dict
            3体間をkey、その3体間のカウントがvaluesとなるdict
        """
        connect_list = self.get_connect_list(
            cut_off=cut_off, bond_type=bond_type)
        # count_triplets[(neibour1, mid, neibour2)] : neibour1-mid-neibour2構造の数
        count_triplets = dict()
        if condition is None:
            target_atoms = np.array([True] * self.get_total_atoms())
        else:
            target_atoms = condition(self)
        atom_type_set = self.get_atom_type_set()
        for atom_type1 in atom_type_set:
            for atom_type2 in atom_type_set:
                for atom_type3 in atom_type_set:
                    count_triplets[(atom_type1, atom_type2, atom_type3)] = 0
        sdat_atom_type = self.atoms['type'].values
        for mid_atom_idx, c_list in enumerate(connect_list):
            if not target_atoms[mid_atom_idx]:
                continue
            mid_atom_type = sdat_atom_type[mid_atom_idx]
            for neibour_atom_idx1, neibour_atom_idx2 in itertools.combinations(c_list, 2):
                if not target_atoms[neibour_atom_idx1]:
                    continue
                if not target_atoms[neibour_atom_idx2]:
                    continue
                neibour_atom_type1 = sdat_atom_type[neibour_atom_idx1]
                neibour_atom_type2 = sdat_atom_type[neibour_atom_idx2]
                if neibour_atom_type1 == neibour_atom_type2:
                    count_triplets[(neibour_atom_type1,
                                    mid_atom_type, neibour_atom_type2)] += 1
                else:
                    count_triplets[(neibour_atom_type1,
                                    mid_atom_type, neibour_atom_type2)] += 1
                    count_triplets[(neibour_atom_type2,
                                    mid_atom_type, neibour_atom_type1)] += 1

        count_triplets_renamed = dict()
        for (neibour1, mid, neibour2), count in count_triplets.items():
            neibour_type1 = self.atom_type_to_symbol[neibour1]
            mid_type = self.atom_type_to_symbol[mid]
            neibour_type2 = self.atom_type_to_symbol[neibour2]
            triplets_renamed = f'{neibour_type1}-{mid_type}-{neibour_type2}'
            count_triplets_renamed[triplets_renamed] = count
        return count_triplets_renamed

    def count_atom_types(self, res_type='series', condition=None):
        """原子のタイプごとに原子の個数をカウントする関数
        Parameters
        ----------
        res_type : str
            res_type='series'のときは結果をpd.Seriesで返す
            res_type='dict'のときは結果をdictで返す
        condition : function
            condition関数
        """
        if condition is None:
            target_atoms = np.array([True] * self.get_total_atoms())
        else:
            target_atoms = condition(self)
        if res_type == 'series':
            return self.atoms.loc[target_atoms,'type'].value_counts().rename(index=self.atom_type_to_symbol)
        elif res_type == 'dict':
            return self.atoms.loc[target_atoms,'type'].value_counts().rename(index=self.atom_type_to_symbol).to_dict()
        else:
            raise ValueError(
                f'res_type: {res_type} is not supported. supported res_type : [series, dict]')

    def get_coordination_number_distribution(self, cut_off=0.5, target_atom_type=-1, bond_type='dumpbond', condition=None) -> dict:
        """
        配位数の分布を調べる
        Parameters
        ----------
        cut_off : float
            bond_type == 'dumppos' の時はcut_offの単位はÅ
            ある原子からcut_off以下の距離にある原子は結合しているとみなす

            bond_type == 'dumpbond' の時はcut_offの単位はbond order
            ある原子とある原子のbond orderの和がcut_off以上のときに結合しているとみなす

            bond_type == 'dumpbond_cg' の時はcut_offは不必要
        target_atom_type : int
            対象となる原子のtype
        bond_type : str
            bond_type == 'dumppos' の時はconnect_listはdumpposから生成される
            bond_type == 'dumpbond' の時connect_listはdumpbondから生成される
            bond_type == 'dumpbond_cg' の時connect_listはdumpbond_cgから生成される

        condition : function
            配位数を調べたい原子を指定する関数

        Returns
        -------
            coordination_number_distribution: dict
                配位数の分布
        """
        assert target_atom_type in self.get_atom_type_set(), 'invalid target_atom_type.'
        if condition is None:
            target_atoms = np.array([True] * self.get_total_atoms())
        else:
            target_atoms = condition(self)
        connect_list = self.get_connect_list(
            cut_off=cut_off, bond_type=bond_type)
        coordination_counter = {i: 0 for i in range(9)}
        atom_types = self.atoms['type'].values
        for atom_idx in range(self.get_total_atoms()):
            if not target_atoms[atom_idx]:
                continue
            if atom_types[atom_idx] != target_atom_type:
                continue
            coordination_counter[len(connect_list[atom_idx])] += 1
        return coordination_counter


class AnalyzeAtomForSDats():
    def __init__():
        pass

    def count_triplets(self, cut_off=0.5, bond_type='dumpbond', condition=None) -> pd.DataFrame:
        """
        ３体間の原子をカウントする関数
        例えば、水分子一つにH-O-Hは1個含まれている。
        メタン分子一つにはcombination(4,2) = 6より、H-C-Hは6個含まれている。

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

        condition : function
            3体間のうち、カウントしたい部分を指定する関数

        Returns
        -------
        df_count_triplets : DataFrame
            step_numをindex、3体間をcolumns、その3体間のカウントがdataとなるDataFrame
        """
        connect_lists = self.get_connect_lists(
            cut_off=cut_off, bond_type=bond_type)
        # count_triplets[step_idx][(neibour1, mid, neibour2)] :step_idxでのneibour1-mid-neibour2構造の数
        count_triplets = [dict() for _ in range(len(self.step_nums))]
        atom_type_set = self.get_atom_type_set()
        for step_idx in range(len(self.step_nums)):
            for atom_type1 in atom_type_set:
                for atom_type2 in atom_type_set:
                    for atom_type3 in atom_type_set:
                        count_triplets[step_idx][(
                            atom_type1, atom_type2, atom_type3)] = 0

        for step_idx in trange(len(self.step_nums), desc='[counting triplets]'):
            if condition is None:
                target_atoms = np.array([True] * self.get_total_atoms())
            else:
                target_atoms = condition(self, step_idx)

            sdat_atom_type = self.atoms[step_idx]['type'].values
            for mid_atom_idx, c_list in enumerate(connect_lists[step_idx]):
                if not target_atoms[mid_atom_idx]:
                    continue
                mid_atom_type = sdat_atom_type[mid_atom_idx]
                for neibour_atom_idx1, neibour_atom_idx2 in itertools.combinations(c_list, 2):
                    if not target_atoms[neibour_atom_idx1]:
                        continue
                    if not target_atoms[neibour_atom_idx2]:
                        continue
                    neibour_atom_type1 = sdat_atom_type[neibour_atom_idx1]
                    neibour_atom_type2 = sdat_atom_type[neibour_atom_idx2]
                    if neibour_atom_type1 == neibour_atom_type2:
                        count_triplets[step_idx][(neibour_atom_type1,
                                                  mid_atom_type, neibour_atom_type2)] += 1
                    else:
                        count_triplets[step_idx][(neibour_atom_type1,
                                                  mid_atom_type, neibour_atom_type2)] += 1
                        count_triplets[step_idx][(neibour_atom_type2,
                                                  mid_atom_type, neibour_atom_type1)] += 1
        df_count_triplets = pd.DataFrame(count_triplets, index=self.step_nums)

        def change_column_name(column):
            atom_type1, atom_type2, atom_type3 = column
            return f'{self.atom_type_to_symbol[atom_type1]}-{self.atom_type_to_symbol[atom_type2]}-{self.atom_type_to_symbol[atom_type3]}'
        df_count_triplets = df_count_triplets.rename(
            columns=change_column_name)

        return df_count_triplets

    def count_atom_types(self):
        count_atom_types_list = []
        for step_idx in range(len(self.step_nums)):
            count_dict = self.atoms[step_idx]['type'].value_counts().to_dict()
            count_atom_types_list.append(count_dict)
        df_count_atom_types = pd.DataFrame(
            data=count_atom_types_list, index=self.step_nums)
        df_count_atom_types.rename(
            columns=self.atom_type_to_symbol, inplace=True)
        return df_count_atom_types

    def get_coordination_number_distribution(self, cut_off=0.5, target_atom_type=-1, bond_type='dumpbond', condition=None) -> pd.DataFrame:
        """
        配位数の分布を調べる
        Parameters
        ----------
        cut_off : float
            bond_type == 'dumppos' の時はcut_offの単位はÅ
            ある原子からcut_off以下の距離にある原子は結合しているとみなす

            bond_type == 'dumpbond' の時はcut_offの単位はbond order
            ある原子とある原子のbond orderの和がcut_off以上のときに結合しているとみなす

            bond_type == 'dumpbond_cg' の時はcut_offは不必要
        target_atom_type : int
            対象となる原子のtype
        bond_type : str
            bond_type == 'dumppos' の時はconnect_listはdumpposから生成される
            bond_type == 'dumpbond' の時connect_listはdumpbondから生成される
            bond_type == 'dumpbond_cg' の時connect_listはdumpbond_cgから生成される

        condition : function
            配位数を調べたい原子を指定する関数

        Returns
        -------
            df_coordination_number_distribution: pd.DataFrame
                配位数の分布
        """
        assert target_atom_type in self.get_atom_type_set(), 'invalid target_atom_type.'
        connect_lists = self.get_connect_lists(
            cut_off=cut_off, bond_type=bond_type)
        coordination_list = []
        for step_idx in trange(len(self.step_nums), desc='[counting coordination num]'):
            if condition is None:
                target_atoms = np.array([True] * self.get_total_atoms())
            else:
                target_atoms = condition(self, step_idx)

            coordination_counter = {i: 0 for i in range(9)}
            atom_types = self.atoms[step_idx]['type'].values
            for atom_idx in range(self.get_total_atoms()):
                if not target_atoms[atom_idx]:
                    continue
                if atom_types[atom_idx] != target_atom_type:
                    continue
                coordination_counter[len(connect_lists[step_idx][atom_idx])] += 1
            coordination_list.append(coordination_counter)
        
        df_coordination_number = pd.DataFrame(data=coordination_list, index=self.step_nums)
        return df_coordination_number