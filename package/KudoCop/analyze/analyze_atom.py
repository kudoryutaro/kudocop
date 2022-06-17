from tqdm import tqdm, trange
import itertools
import pandas as pd
import numpy as np


class AnalyzeAtom():
    def __init__():
        pass

    def count_triplets(self, cut_off=0.5,bond_type='dumpbond', condition=None):
        connect_list = self.get_connect_list(cut_off=cut_off,bond_type=bond_type)
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
                count_triplets[(neibour_atom_type1,
                                mid_atom_type, neibour_atom_type2)] += 1
        count_triplets_renamed = dict()
        for (neibour1, mid, neibour2), count in count_triplets.items():
            neibour_type1 = self.atom_type_to_symbol[neibour1]
            mid_type = self.atom_type_to_symbol[mid]
            neibour_type2 = self.atom_type_to_symbol[neibour2]
            triplets_renamed = f'{neibour_type1}-{mid_type}-{neibour_type2}'
            count_triplets_renamed[triplets_renamed] = count
        return count_triplets_renamed


class AnalyzeAtomForSDats():
    def __init__():
        pass

    def count_triplets(self, cut_off: float, condition=None):
        self.get_connect_lists(cut_off)
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
            for mid_atom_idx, c_list in enumerate(self.connect_lists[step_idx]):
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
                    count_triplets[step_idx][(neibour_atom_type1,
                                              mid_atom_type, neibour_atom_type2)] += 1

        df_count_triplets = pd.DataFrame(count_triplets, index=self.step_nums)

        def change_column_name(column):
            atom_type1, atom_type2, atom_type3 = column
            return f'{self.atom_type_to_symbol[atom_type1]}-{self.atom_type_to_symbol[atom_type2]}-{self.atom_type_to_symbol[atom_type3]}'
        df_count_triplets = df_count_triplets.rename(
            columns=change_column_name)

        return df_count_triplets
