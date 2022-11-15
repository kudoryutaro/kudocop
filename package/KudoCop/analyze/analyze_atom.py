from tqdm import tqdm, trange
import itertools
import pandas as pd
import numpy as np

try:
    from . import neighbor
    import ase
    from ase.atoms import Atoms
    import torch
except:
    pass

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


    def count_structure_target_group(self, cut_off=0.5, target_atom_type=-1, bond_type='dumpbond', condition=None) -> dict:
        """
        ある原子種をtargetとして、targetに結合している原子の種類と数のタプルを返す
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
            count_structure_target_group: dict
                targetに結合している原子種と数
        """
        count_structure_target_group = {}
        connect_list = self.get_connect_list(bond_type=bond_type, cut_off=cut_off)
        type_max_num = max(self.get_atom_type_set())
        atom_types = self.atoms['type'].values
        for atom_idx in trange(self.get_total_atoms()):
            if atom_types[atom_idx] != target_atom_type:
                continue
            current_neighbor_atom_type_count = [0] * (type_max_num + 1)
            for neighbor_atom_idx in connect_list[atom_idx]:
                neighbor_atom_type = atom_types[neighbor_atom_idx]
                current_neighbor_atom_type_count[neighbor_atom_type] += 1
            current_neighbor_atom_type_count = tuple(current_neighbor_atom_type_count)
            if current_neighbor_atom_type_count not in count_structure_target_group:
                count_structure_target_group[current_neighbor_atom_type_count] = 0
            count_structure_target_group[current_neighbor_atom_type_count] += 1
        return count_structure_target_group

    # def get_relative_coordinates_info(self, cut_off:float):
    #     """それぞれの原子から見て、距離がcut_off以下の原子を探し、相対座標, 距離, タイプ, atom_idxを返す関数
    #     Parameters
    #     ----------
    #         cut_off : float
    #             カットオフ 単位はÅ

    #     Returns
    #     -------
    #         connect_list : list
    #             連結リスト
    #         relative_coordinates : list
    #             relative_coordinates[atom_idx] = [[x_ij, y_ij, z_ij],
    #                                               [x_ij, y_ij, z_ij],
    #                                               [x_ij, y_ij, z_ij],
    #                                               [x_ij, y_ij, z_ij]]
    #     Note
    #     ----
    #         connect_list[atom_idx]中のインデックスと
    #         relative_coordinates[atom_idx]中のインデックスは一致している
    #         Example
    #         -------
    #             connect_list[3] = [0, 1, 4]
    #             relative_coordinates[atom_idx] = [[2.1, 4.2, 2.3],
    #                                               [6.3, -3.1, 1.8],
    #                                               [2.1, 1.2, -2.1]]
    #             の場合はatomのインデックスが3の原子の周りにatomのインデックスが0, 1, 4の原子がある
    #             atomのインデックスが3の原子から見て、atomのインデックスが0の原子は[2.1, 4.2, 2.3]方向にある
    #             atomのインデックスが3の原子から見て、atomのインデックスが1の原子は[6.3, -3.1, 1.8]方向にある
    #             atomのインデックスが3の原子から見て、atomのインデックスが4の原子は[2.1, 1.2, -2.1]方向にある

    #     """
    #     ase_atoms = Atoms(
    #         positions=self.atoms[['x','y','z']],
    #         symbols=self.atoms['type'].map(self.atom_type_to_symbol),
    #         cell=self.cell,
    #         pbc=[1, 1, 1])
    #     atom_types = self.atoms['type'].values
    #     ovito_data = ase_to_ovito(ase_atoms)
    #     ovito_neighbor_finder = ovito.data.CutoffNeighborFinder(cut_off, ovito_data)
    #     relative_coordinates_list = [[] for _ in range(self.get_total_atoms())]
    #     relative_atom_type_list = [[] for _ in range(self.get_total_atoms())]
    #     relative_dist_list = [[] for _ in range(self.get_total_atoms())]
    #     connect_list = [[] for _ in range(self.get_total_atoms())]
    #     for atom_idx in range(self.get_total_atoms()):
    #         for neighbor_atom in ovito_neighbor_finder.find(atom_idx):
    #             relative_coordinates_list[atom_idx].append(neighbor_atom.index)
    #             relative_dist_list[atom_idx].append(neighbor_atom.distance)
    #             relative_coordinates_list[atom_idx].append(neighbor_atom.delta)
    #             relative_atom_type_list[atom_idx].append(atom_types[neighbor_atom.index])
    #             connect_list[atom_idx].append(neighbor_atom.index)

    #     return connect_list, relative_coordinates_list, relative_atom_type_list, relative_dist_list

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

    def count_atom_types(self, condition=None):
        count_atom_types_list = []
        for step_idx in range(len(self.step_nums)):
            if condition is None:
                target_atoms = np.array([True] * self.get_total_atoms())
            else:
                target_atoms = condition(self, step_idx)
            count_dict = self.atoms[step_idx].loc[target_atoms, 'type'].value_counts().to_dict()
            count_atom_types_list.append(count_dict)
        df_count_atom_types = pd.DataFrame(
            data=count_atom_types_list, index=self.step_nums).fillna(0).astype(int)
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

            coordination_counter = {i: 0 for i in range(13)}
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

    
    def get_time_variation_of_coordinated_atoms(self, cut_off=0.5, bond_type='dumpbond', 
                                target_atom_type=-1, condition=None, rename=True):
        """ある種類の原子に注目して、その原子に何が配位しているかを調べる。そして、その原子が次のステップでどのような
        配位に変わったのかをカウントする。
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
        rename : bool
            rename=Trueにすると読みやすい形式で返す
        Returns
        -------
            time_variation_of_coordinated_atoms: dict
                配位数の時間変化
        Explanation
        -----------
            time_variation_of_coordinated_atomsの形式
                time_variation_of_coordinated_atoms[今のステップの構造][次のステップの構造] = その構造のカウント数
        Example
        -------
            0, 1000, 2000ステップ読み込んでいるとし、target_atom_type=4 (N, 窒素)とする
            この時、0->1000ステップ、1000->2000ステップのNの配位数の変化がカウントされる
            >>>sdats.atom_type_to_symbol
            {1: 'C', 2: 'H', 3: 'O', 4: 'N', 5: 'S', 6: 'Si', 7: 'Na', 8: 'F', 9: 'P'}
            >>>time_variation_of_coordinated_atoms = sdats.get_time_variation_of_coordinated_atoms(target_atom_type=4)
            この例では、H2つとSi1つと結合しているNは次のステップで、H1つ, O1つ, Si1つと結合するようになるNは62個あった。
            つまり、'N-H2Si1' -> 'N-H1O1Si1' という変化をしたNの個数は62個あった。
            >>>pprint(time_variation_of_coordinated_atoms['N-H2Si1'])
            {'N-H1O1Si1': 62,
            'N-H1Si1': 67,
            'N-H1Si2': 177,
            'N-H2O1Si1': 70,
            'N-H2Si1': 1663,
            'N-H2Si2': 145,
            'N-H3': 184,
            'N-H3Si1': 277}
        """
        assert target_atom_type in self.get_atom_type_set(), 'invalid target_atom_type.'
        assert len(self.step_nums) >= 2, 'ステップ数が少なすぎます'
        connect_lists = self.get_connect_lists(
            cut_off=cut_off, bond_type=bond_type)
        atom_types = self.atoms[0]['type'].values
        atom_type_max_num = max(atom_types)
        # time_variation_of_coordinated_atoms[coordination count of current_step][coordination count of next_step] = count 
        time_variation_of_coordinated_atoms = dict()

        for current_step_idx in trange(len(self.step_nums) - 1, desc='[counting time var of coord atoms]'):
            if condition is None:
                target_atoms = np.array([True] * self.get_total_atoms())
            else:
                target_atoms = condition(self, current_step_idx)
            next_step_idx = current_step_idx + 1
            for atom_idx in range(self.get_total_atoms()):
                if not target_atoms[atom_idx]:
                    continue
                if atom_types[atom_idx] != target_atom_type:
                    continue
                current_coord_counter = [0] * (atom_type_max_num + 1)
                for neighbor_atom_idx in connect_lists[current_step_idx][atom_idx]:
                    neighbor_atom_type = atom_types[neighbor_atom_idx]
                    current_coord_counter[neighbor_atom_type] += 1
                current_coord_counter = tuple(current_coord_counter)

                next_coord_counter = [0] * (atom_type_max_num + 1)
                for neighbor_atom_idx in connect_lists[next_step_idx][atom_idx]:
                    neighbor_atom_type = atom_types[neighbor_atom_idx]
                    next_coord_counter[neighbor_atom_type] += 1
                next_coord_counter = tuple(next_coord_counter)

                if not current_coord_counter in time_variation_of_coordinated_atoms:
                    time_variation_of_coordinated_atoms[current_coord_counter] = {}
                if not next_coord_counter in time_variation_of_coordinated_atoms[current_coord_counter]:
                    time_variation_of_coordinated_atoms[current_coord_counter][next_coord_counter] = 0
                time_variation_of_coordinated_atoms[current_coord_counter][next_coord_counter] += 1
        if not rename:
            return time_variation_of_coordinated_atoms        

        time_variation_of_coordinated_atoms_renamed = dict()
        for current_coord_counter, current_coord_counter_counter in time_variation_of_coordinated_atoms.items():
            current_coord_counter_renamed_list = []
            current_coord_counter_renamed_list.append(f'{self.atom_type_to_symbol[target_atom_type]}-')
            for atom_type in range(1, atom_type_max_num + 1):
                if current_coord_counter[atom_type] == 0:
                    continue
                current_coord_counter_renamed_list.append(f'{self.atom_type_to_symbol[atom_type]}{current_coord_counter[atom_type]}')
            current_coord_counter_renamed = ''.join(current_coord_counter_renamed_list)
            if current_coord_counter_renamed not in time_variation_of_coordinated_atoms_renamed:
                time_variation_of_coordinated_atoms_renamed[current_coord_counter_renamed] = {}

            for next_coord_counter, count_num in current_coord_counter_counter.items():
                next_coord_counter_renamed_list = []
                next_coord_counter_renamed_list.append(f'{self.atom_type_to_symbol[target_atom_type]}-')
                for atom_type in range(1, atom_type_max_num + 1):
                    if next_coord_counter[atom_type] == 0:
                        continue
                    next_coord_counter_renamed_list.append(f'{self.atom_type_to_symbol[atom_type]}{next_coord_counter[atom_type]}')
                next_coord_counter_renamed = ''.join(next_coord_counter_renamed_list)
                time_variation_of_coordinated_atoms_renamed[current_coord_counter_renamed]

                time_variation_of_coordinated_atoms_renamed[current_coord_counter_renamed][next_coord_counter_renamed] = count_num

        return time_variation_of_coordinated_atoms_renamed
